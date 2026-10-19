# Internal EM machinery for fitHsmmEM().
#
# The expanded HMM is never built as a dense matrix. Its transition matrix has
# only two kinds of nonzero entries,
#
#   stay:  (i, r) -> (i, r + 1)   [or (i, m_i) -> (i, m_i)]   weight 1 - c_i(r)
#   exit:  (i, r) -> (j, 1)                                    weight omega(i, j) c_i(r)
#
# so a chain is fully described by the concatenated hazard vector, the
# aggregate start/end indices, and omega. Each forward/backward step then costs
# O(N + K^2) instead of the O(N^2) dense product used by scaledForward().
#
# Hazards are always obtained from makeDwellBlock(), i.e. through exactly the
# same code path hsmmToHmm() uses, so the likelihood EM climbs is by
# construction the one fitHsmm() optimizes with nlm().


# ---------------------------------------------------------------------------
# Chain structure
# ---------------------------------------------------------------------------

buildChainStructure <- function(model) {
  nStates <- model$nStates

  hazardList <- lapply(seq_len(nStates), function(stateIndex) {
    makeDwellBlock(
      sojournParameters = model$sojournParameters[[stateIndex]],
      sojournType = model$sojournType[stateIndex]
    )$hazard
  })

  sizes <- lengths(hazardList)
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  nSub <- sum(sizes)
  aggregate <- rep.int(seq_len(nStates), sizes)
  hazard <- unlist(hazardList, use.names = FALSE)

  # Successor under "stay": the next substate, or the last substate itself.
  stayTo <- seq_len(nSub) + 1L
  stayTo[ends] <- ends

  # N x K indicator matrix: collapses substate columns to HSMM-state columns.
  aggregateIndicator <- matrix(0, nrow = nSub, ncol = nStates)
  aggregateIndicator[cbind(seq_len(nSub), aggregate)] <- 1

  list(
    nStates = nStates,
    nSub = nSub,
    sizes = sizes,
    starts = starts,
    ends = ends,
    aggregate = aggregate,
    hazard = hazard,
    stay = 1 - hazard,
    stayTo = stayTo,
    aggregateIndicator = aggregateIndicator,
    omega = model$tpm,
    hazardList = hazardList
  )
}


# Stationary distribution of the expanded chain, in closed form.
#
# If e_i is the long-run rate of entering state i, then
#   pi(i, r)   = e_i * S_i(r),            r < m_i,  S_i(r) = prod_{s < r} (1 - c_i(s))
#   pi(i, m_i) = e_i * S_i(m_i) / c_i(m_i)  (the self-loop accumulates the tail)
# and the exit flux from i equals e_i, so e is proportional to the stationary
# distribution of omega. O(N + K^3) instead of an N x N linear solve.
expandedStationaryDistribution <- function(chain) {
  entryRates <- getStationaryDistribution(chain$omega)

  pieces <- lapply(seq_len(chain$nStates), function(stateIndex) {
    hazard <- chain$hazardList[[stateIndex]]
    m <- length(hazard)
    survival <- c(1, cumprod(1 - hazard))[seq_len(m)]

    if (hazard[m] <= 0) {
      stop("Tail hazard of state ", stateIndex, " is zero; no stationary distribution.")
    }

    survival[m] <- survival[m] / hazard[m]
    entryRates[stateIndex] * survival
  })

  delta <- unlist(pieces, use.names = FALSE)
  delta / sum(delta)
}


getExpandedInitial <- function(model, chain, stationary) {
  if (stationary) {
    return(expandedStationaryDistribution(chain))
  }

  if (is.null(model$delta)) {
    stop("model$delta is required when stationary = FALSE.")
  }

  delta <- numeric(chain$nSub)
  delta[chain$starts] <- model$delta / sum(model$delta)
  delta
}


# ---------------------------------------------------------------------------
# Sparse scaled forward-backward
# ---------------------------------------------------------------------------

# One step of alpha %*% Gamma without forming Gamma.
propagateForward <- function(alphaRow, chain) {
  nSub <- chain$nSub
  starts <- chain$starts
  ends <- chain$ends

  # Stay: mass at (i, r) moves to (i, r + 1); the last substate also self-loops.
  stayMass <- alphaRow * chain$stay
  predicted <- c(0, stayMass[-nSub])
  predicted[starts] <- 0
  predicted[ends] <- predicted[ends] + stayMass[ends]

  # Exit: F(i) = sum_r alpha(i, r) c_i(r), then redistribute through omega.
  cumulativeExit <- cumsum(alphaRow * chain$hazard)
  exitFlux <- diff(c(0, cumulativeExit[ends]))
  predicted[starts] <- predicted[starts] + as.vector(exitFlux %*% chain$omega)

  predicted
}


# Scaled forward-backward on the expanded chain.
#
# Emissions enter as exp(logB - rowMax), i.e. each time step's densities are
# divided by their largest value. The recursion is invariant to a per-time
# constant, so this prevents the whole-row underflow that raw densities suffer
# for very negative log-likelihoods (e.g. curve-valued observations), and the
# offsets are added back to the log-likelihood exactly.
#
# Conventions: alphaHat rows sum to one; gamma = alphaHat * betaHat;
# logLik = sum(log(scale)) + sum(offset).
forwardBackwardExpanded <- function(chain, delta, logEmission, relativeFloor = 1e-300) {
  nObs <- nrow(logEmission)
  nSub <- chain$nSub
  aggregate <- chain$aggregate

  offset <- apply(logEmission, 1, max)

  if (any(!is.finite(offset))) {
    return(NULL)
  }

  scaledEmission <- pmax(exp(logEmission - offset), relativeFloor)

  alphaHat <- matrix(0, nrow = nObs, ncol = nSub)
  betaHat <- matrix(0, nrow = nObs, ncol = nSub)
  scale <- numeric(nObs)

  raw <- delta * scaledEmission[1, aggregate]
  scale[1] <- sum(raw)
  if (!is.finite(scale[1]) || scale[1] <= 0) return(NULL)
  alphaHat[1, ] <- raw / scale[1]

  if (nObs >= 2L) {
    for (timeIndex in 2:nObs) {
      raw <- propagateForward(alphaHat[timeIndex - 1L, ], chain) *
        scaledEmission[timeIndex, aggregate]
      scale[timeIndex] <- sum(raw)
      if (!is.finite(scale[timeIndex]) || scale[timeIndex] <= 0) return(NULL)
      alphaHat[timeIndex, ] <- raw / scale[timeIndex]
    }
  }

  betaHat[nObs, ] <- 1

  if (nObs >= 2L) {
    for (timeIndex in seq(nObs - 1L, 1L)) {
      # u(k) = b_{t+1}(k) * betaHat_{t+1}(k)
      nextValue <- scaledEmission[timeIndex + 1L, aggregate] * betaHat[timeIndex + 1L, ]
      # G(i) = sum_{j != i} omega(i, j) u(j, 1)
      exitValue <- as.vector(chain$omega %*% nextValue[chain$starts])
      betaHat[timeIndex, ] <- (
        chain$stay * nextValue[chain$stayTo] + chain$hazard * exitValue[aggregate]
      ) / scale[timeIndex + 1L]
    }
  }

  list(
    alphaHat = alphaHat,
    betaHat = betaHat,
    scale = scale,
    scaledEmission = scaledEmission,
    logLik = sum(log(scale)) + sum(offset)
  )
}


# ---------------------------------------------------------------------------
# E-step
# ---------------------------------------------------------------------------

#' E-step for the expanded-HMM EM algorithm (internal)
#'
#' Runs the sparse forward-backward pass and returns the expected sufficient
#' statistics used by the M-step:
#' \describe{
#'   \item{`atRisk`}{eta_i(r) = sum_{t < T} gamma_t(i, r), per substate}
#'   \item{`departures`}{d_i(r) = expected exits from (i, r), per substate}
#'   \item{`transitions`}{K x K matrix N(i, j) of expected i -> j jumps}
#'   \item{`stateWeights`}{T x K matrix gamma_t(i) = sum_r gamma_t(i, r)}
#'   \item{`gammaFirst`}{gamma_1(i, r), for the non-stationary delta update}
#' }
#' All statistics are formed with a handful of matrix products after the
#' passes; no T x N x N xi array is ever stored. `prepared` is the emission
#' family's reduced data (see emissionFamilies.R); pass it to avoid
#' recomputing it every iteration.
#' @noRd
hsmmEStep <- function(model, x, stationary, prepared = NULL) {
  chain <- buildChainStructure(model)
  delta <- getExpandedInitial(model, chain, stationary)

  family <- getEmissionFamily(model$emissionType)
  if (is.null(prepared)) {
    prepared <- family$prepare(x, model$emissionControl)
  }
  logEmission <- family$logDensity(prepared, model$emissionParameters, model$emissionControl)

  fb <- forwardBackwardExpanded(chain, delta, logEmission)

  if (is.null(fb)) {
    return(list(logLik = -Inf, chain = chain))
  }

  nObs <- nrow(logEmission)
  notLast <- seq_len(nObs - 1L)
  scaleNext <- fb$scale[-1L]

  gamma <- fb$alphaHat * fb$betaHat                        # T x N, rows sum to 1
  nextValue <- fb$scaledEmission[, chain$aggregate] * fb$betaHat
  entryValue <- nextValue[-1L, chain$starts, drop = FALSE] / scaleNext   # (T-1) x K

  # G_t(i) / c_{t+1}, one row per t = 1, ..., T - 1.
  exitValue <- entryValue %*% t(chain$omega)

  # d_i(r) = c_i(r) * sum_t alphaHat_t(i, r) G_t(i) / c_{t+1}
  departures <- chain$hazard * colSums(
    fb$alphaHat[notLast, , drop = FALSE] * exitValue[, chain$aggregate, drop = FALSE]
  )

  # F_t(i) = sum_r alphaHat_t(i, r) c_i(r);  N(i, j) = omega(i, j) sum_t F_t(i) u_{t+1}(j, 1) / c_{t+1}
  exitFlux <- sweep(fb$alphaHat[notLast, , drop = FALSE], 2, chain$hazard, `*`) %*%
    chain$aggregateIndicator
  transitions <- chain$omega * crossprod(exitFlux, entryValue)

  list(
    logLik = fb$logLik,
    chain = chain,
    atRisk = colSums(gamma[notLast, , drop = FALSE]),
    departures = departures,
    transitions = transitions,
    stateWeights = gamma %*% chain$aggregateIndicator,
    gammaFirst = gamma[1, ]
  )
}


# ---------------------------------------------------------------------------
# M-step: dwell times
# ---------------------------------------------------------------------------

# x * log(y) with the convention 0 * log(0) = 0.
xLogY <- function(x, y) {
  out <- numeric(length(x))
  positive <- x > 0
  out[positive] <- x[positive] * log(y[positive])
  out
}


# Q_dwell for one state: weighted Bernoulli log-likelihood of the hazards.
dwellObjective <- function(hazard, departures, atRisk) {
  sum(xLogY(departures, hazard) + xLogY(atRisk - departures, 1 - hazard))
}


hazardToPmf <- function(hazard) {
  survival <- c(1, cumprod(1 - hazard))[seq_along(hazard)]
  hazard * survival
}


clampProb <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}


# Closed-form life-table hazard estimate, keeping the old hazard wherever the
# risk set is empty.
lifeTableHazard <- function(departures, atRisk, hazardOld, minCount) {
  hazard <- hazardOld
  informative <- atRisk > minCount
  hazard[informative] <- pmin(pmax(departures[informative] / atRisk[informative], 0), 1)
  hazard
}


# Fisher scoring on the working scale for any parametric family.
#
#   score U = sum_r (d_r - eta_r h_r) / (h_r (1 - h_r)) * dh_r
#   info  I = sum_r eta_r / (h_r (1 - h_r)) * dh_r dh_r'
#
# dh is taken by central differences through makeDwellBlock(), so the same
# routine serves pois, shiftp, nb, shiftnb and weibull, and future families
# need no new derivative code. Step-halving guarantees Q_dwell never decreases
# (generalized EM), so a few steps per EM iteration are enough.
fisherScoringDwell <- function(
    sojournParameters,
    sojournType,
    departures,
    atRisk,
    maxSteps = 5,
    stepTol = 1e-8,
    finiteDiffStep = 1e-5,
    maxStepSize = 3,
    ridge = 1e-10
) {
  hazardAt <- function(working) {
    constrained <- workingToSojournOne(
      workingVec = working,
      sojournTemplate = sojournParameters,
      sojournType = sojournType,
      fitShift = FALSE
    )
    makeDwellBlock(constrained, sojournType)$hazard
  }

  working <- sojournToWorkingOne(sojournParameters, sojournType, fitShift = FALSE)
  hazard <- hazardAt(working)
  objective <- dwellObjective(hazard, departures, atRisk)
  nPar <- length(working)

  for (stepIndex in seq_len(maxSteps)) {
    jacobian <- matrix(0, nrow = length(hazard), ncol = nPar)

    for (parIndex in seq_len(nPar)) {
      step <- finiteDiffStep * max(1, abs(working[parIndex]))
      up <- working
      down <- working
      up[parIndex] <- up[parIndex] + step
      down[parIndex] <- down[parIndex] - step
      jacobian[, parIndex] <- (hazardAt(up) - hazardAt(down)) / (2 * step)
    }

    # Durations with h = 0 (below a shift) or h = 1, or with no risk set, carry
    # no information and would give 0/0 weights.
    valid <- hazard > 1e-12 & hazard < 1 - 1e-12 & atRisk > 0

    if (!any(valid)) {
      break
    }

    variance <- hazard[valid] * (1 - hazard[valid])
    residual <- (departures[valid] - atRisk[valid] * hazard[valid]) / variance
    score <- colSums(residual * jacobian[valid, , drop = FALSE])
    information <- crossprod(jacobian[valid, , drop = FALSE] * sqrt(atRisk[valid] / variance))

    direction <- tryCatch(
      solve(information + ridge * diag(nPar), score),
      error = function(e) NULL
    )

    if (is.null(direction) || any(!is.finite(direction))) {
      break
    }

    largest <- max(abs(direction))
    if (largest > maxStepSize) {
      direction <- direction * (maxStepSize / largest)
    }

    accepted <- FALSE

    for (halving in 0:30) {
      candidate <- working + direction / 2^halving
      candidateHazard <- tryCatch(hazardAt(candidate), error = function(e) NULL)

      if (is.null(candidateHazard)) next

      candidateObjective <- dwellObjective(candidateHazard, departures, atRisk)

      if (is.finite(candidateObjective) && candidateObjective >= objective) {
        accepted <- TRUE
        break
      }
    }

    if (!accepted) {
      break
    }

    moved <- max(abs(candidate - working))
    working <- candidate
    hazard <- candidateHazard
    objective <- candidateObjective

    if (moved < stepTol) {
      break
    }
  }

  workingToSojournOne(
    workingVec = working,
    sojournTemplate = sojournParameters,
    sojournType = sojournType,
    fitShift = FALSE
  )
}


updateDwellOne <- function(
    sojournParameters,
    sojournType,
    hazardOld,
    departures,
    atRisk,
    control
) {
  m <- length(hazardOld)

  if (sojournType == "np") {
    hazard <- lifeTableHazard(departures, atRisk, hazardOld, control$minCount)
    hazard[m] <- 1                 # finite support
    return(hazardToPmf(hazard))
  }

  if (sojournType == "npg") {
    hazard <- lifeTableHazard(departures, atRisk, hazardOld, control$minCount)
    hazard <- clampProb(hazard)    # tail hazard in (0, 1)
    headPmf <- hazardToPmf(hazard)[-m]

    # Several head hazards near one (every observed sojourn ended by some
    # r < m) make 1 - prod(1 - h) round to exactly 1, which npg cannot store.
    maxHeadMass <- 1 - 1e-8
    if (sum(headPmf) > maxHeadMass) {
      headPmf <- headPmf * (maxHeadMass / sum(headPmf))
    }

    return(c(headPmf, hazard[m]))  # npg stores the tail *hazard* last
  }

  if (sojournType == "geo") {
    totalRisk <- sum(atRisk)
    if (totalRisk <= control$minCount) return(sojournParameters)
    prob <- clampProb(sum(departures) / totalRisk)
    names(prob) <- names(sojournParameters)
    return(prob)
  }

  if (sojournType %in% c("pois", "shiftp", "nb", "shiftnb", "weibull")) {
    if (sum(atRisk) <= control$minCount) return(sojournParameters)

    updated <- fisherScoringDwell(
      sojournParameters = sojournParameters,
      sojournType = sojournType,
      departures = departures,
      atRisk = atRisk,
      maxSteps = control$dwellSteps
    )
    names(updated) <- names(sojournParameters)
    return(updated)
  }

  stop("Unsupported sojournType: ", sojournType)
}


# ---------------------------------------------------------------------------
# M-step: transitions, emissions, initial distribution
# ---------------------------------------------------------------------------

updateEmbeddedTpm <- function(tpm, transitions, minCount) {
  nStates <- nrow(tpm)

  # With K = 2 the zero-diagonal constraint fixes omega completely.
  if (nStates <= 2L) {
    return(tpm)
  }

  updated <- tpm

  for (stateIndex in seq_len(nStates)) {
    counts <- transitions[stateIndex, -stateIndex]

    if (sum(counts) > minCount) {
      updated[stateIndex, -stateIndex] <- counts / sum(counts)
    }
  }

  diag(updated) <- 0
  updated
}


#' M-step for the expanded-HMM EM algorithm (internal)
#'
#' The expected complete-data log-likelihood separates into omega, per-state
#' dwell, per-state emission and initial-distribution blocks, each maximized
#' (or, for parametric dwell families, increased) independently.
#'
#' With `stationary = TRUE` the initial distribution is a function of the
#' other parameters; it is not updated here but recomputed from them in the
#' next E-step (plug-in approximation, so EM is then not exact).
#' @noRd
hsmmMStep <- function(model, eStep, prepared, stationary, control) {
  chain <- eStep$chain
  nStates <- model$nStates

  tpm <- updateEmbeddedTpm(model$tpm, eStep$transitions, control$minCount)

  sojournParameters <- vector("list", nStates)

  for (stateIndex in seq_len(nStates)) {
    index <- chain$starts[stateIndex]:chain$ends[stateIndex]

    sojournParameters[[stateIndex]] <- updateDwellOne(
      sojournParameters = model$sojournParameters[[stateIndex]],
      sojournType = model$sojournType[stateIndex],
      hazardOld = chain$hazard[index],
      departures = eStep$departures[index],
      atRisk = eStep$atRisk[index],
      control = control
    )
  }

  emissionParameters <- getEmissionFamily(model$emissionType)$update(
    prepared,
    eStep$stateWeights,
    model$emissionParameters,
    model$emissionControl,
    control
  )

  delta <- model$delta

  if (!stationary) {
    delta <- eStep$gammaFirst[chain$starts]
    delta <- pmax(delta, 0) / sum(delta)
  }

  buildHsmm(
    tpm = tpm,
    sojournParameters = sojournParameters,
    sojournType = model$sojournType,
    emissionParameters = emissionParameters,
    emissionType = model$emissionType,
    emissionControl = model$emissionControl,
    delta = delta,
    stationary = stationary,
    fitShift = FALSE
  )
}


# ---------------------------------------------------------------------------
# Decoding
# ---------------------------------------------------------------------------

#' Viterbi decoding on the sparse expanded chain, in log space
#'
#' Same recursion as [viterbiLog()] but exploiting the two edge types of the
#' expanded transition matrix, so each step costs `O(N + K^2)`, and taking
#' log-densities directly from the emission family. The latter matters for
#' curve-valued data, whose log-densities are routinely far below
#' `log(.Machine$double.xmin)`: [hsmmViterbi()] exponentiates them first and
#' would underflow.
#'
#' @param model A constrained HSMM model (any emission family).
#' @param x Observations in the format the emission family expects.
#' @param stationary Logical; defaults to `model$stationary`.
#'
#' @return A list with `states` (HSMM states), `expandedStates`, and
#'   `logProb` (log joint probability of the decoded path and the data).
#' @export
hsmmViterbiEM <- function(model, x, stationary = isTRUE(model$stationary)) {
  family <- getEmissionFamily(model$emissionType)
  prepared <- family$prepare(x, model$emissionControl)
  logEmission <- family$logDensity(prepared, model$emissionParameters, model$emissionControl)

  chain <- buildChainStructure(model)
  delta <- getExpandedInitial(model, chain, stationary)

  nObs <- nrow(logEmission)
  nStates <- chain$nStates
  nSub <- chain$nSub
  aggregate <- chain$aggregate
  starts <- chain$starts
  ends <- chain$ends

  logStay <- log(chain$stay)
  logHazard <- log(chain$hazard)
  logOmega <- log(chain$omega)
  advanceFrom <- seq_len(nSub) - 1L        # (i, r - 1) -> (i, r)
  advanceTargets <- setdiff(seq_len(nSub), starts)

  score <- log(delta) + logEmission[1, aggregate]
  backPointer <- matrix(NA_integer_, nrow = nObs, ncol = nSub)

  if (nObs >= 2L) {
    for (timeIndex in 2:nObs) {
      stayScore <- score + logStay
      best <- rep(-Inf, nSub)
      from <- rep(NA_integer_, nSub)

      best[advanceTargets] <- stayScore[advanceFrom[advanceTargets]]
      from[advanceTargets] <- advanceFrom[advanceTargets]

      selfLoop <- stayScore[ends] > best[ends]
      best[ends[selfLoop]] <- stayScore[ends[selfLoop]]
      from[ends[selfLoop]] <- ends[selfLoop]

      exitScore <- score + logHazard
      exitBest <- numeric(nStates)
      exitFrom <- integer(nStates)
      for (i in seq_len(nStates)) {
        segment <- starts[i]:ends[i]
        r <- which.max(exitScore[segment])
        exitBest[i] <- exitScore[segment][r]
        exitFrom[i] <- segment[r]
      }

      entryScore <- exitBest + logOmega            # [i, j]: leave i, enter j
      bestSource <- max.col(t(entryScore), ties.method = "first")
      entryBest <- entryScore[cbind(bestSource, seq_len(nStates))]

      entry <- entryBest > best[starts]
      best[starts[entry]] <- entryBest[entry]
      from[starts[entry]] <- exitFrom[bestSource[entry]]

      score <- best + logEmission[timeIndex, aggregate]
      backPointer[timeIndex, ] <- from
    }
  }

  if (!any(is.finite(score))) {
    stop("Viterbi found no path with positive probability.")
  }

  path <- integer(nObs)
  path[nObs] <- which.max(score)
  if (nObs >= 2L) {
    for (timeIndex in seq(nObs, 2L)) {
      path[timeIndex - 1L] <- backPointer[timeIndex, path[timeIndex]]
    }
  }

  list(states = aggregate[path], expandedStates = path, logProb = max(score))
}


#' Posterior state probabilities
#'
#' @inheritParams hsmmViterbiEM
#' @return An nObs x K matrix of P(S_t = i | data), summed over durations.
#' @export
hsmmStatePosterior <- function(model, x, stationary = isTRUE(model$stationary)) {
  eStep <- hsmmEStep(model, x, stationary)
  if (!is.finite(eStep$logLik)) stop("Model has zero likelihood for these data.")
  eStep$stateWeights
}
