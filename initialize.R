# Internal initialization helpers used by initializeHsmmParameters() and
# initializeHsmmMultiStart() below. None of these are exported.

initializeEmissionParameters <- function(x, nStates, emissionType, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (emissionType == "normal") {
    clusterFit <- kmeans(
      x = matrix(x, ncol = 1),
      centers = nStates,
      nstart = 20
    )
    
    clusterOrder <- order(clusterFit$centers[, 1])
    clusterMap <- match(clusterFit$cluster, clusterOrder)
    
    emissionParameters <- matrix(NA_real_, nrow = nStates, ncol = 2)
    
    for (stateIndex in seq_len(nStates)) {
      stateX <- x[clusterMap == stateIndex]
      
      if (length(stateX) < 2L) {
        stateX <- x
      }
      
      emissionParameters[stateIndex, 1] <- mean(stateX)
      emissionParameters[stateIndex, 2] <- max(sd(stateX), 1e-3)
    }
    
    colnames(emissionParameters) <- c("mean", "sd")
    return(emissionParameters)
  }
  
  if (emissionType == "pois") {
    clusterFit <- kmeans(
      x = matrix(x, ncol = 1),
      centers = nStates,
      nstart = 20
    )
    
    clusterOrder <- order(clusterFit$centers[, 1])
    clusterMap <- match(clusterFit$cluster, clusterOrder)
    
    emissionParameters <- matrix(NA_real_, nrow = nStates, ncol = 1)
    
    for (stateIndex in seq_len(nStates)) {
      stateX <- x[clusterMap == stateIndex]
      
      if (length(stateX) < 1L) {
        stateX <- x
      }
      
      emissionParameters[stateIndex, 1] <- max(mean(stateX), 1e-3)
    }
    
    colnames(emissionParameters) <- "lambda"
    return(emissionParameters)
  }
  
  stop("Unsupported emissionType: ", emissionType)
}


initializeEmbeddedTpm <- function(nStates) {
  if (nStates < 2L) {
    stop("nStates must be at least 2.")
  }
  
  tpm <- matrix(1 / (nStates - 1), nrow = nStates, ncol = nStates)
  diag(tpm) <- 0
  tpm <- tpm / rowSums(tpm)
  
  tpm
}


initializeSojournOne <- function(
    sojournType,
    dwellSize = 10,
    meanDwell = NULL,
    shift = 1,
    fitShift = FALSE
) {
  if (is.null(meanDwell)) {
    meanDwell <- max(2, dwellSize / 2)
  }
  
  if (sojournType == "geo") {
    prob <- min(max(1 / meanDwell, 1e-3), 1 - 1e-3)
    return(c(prob = prob))
  }
  
  if (sojournType == "pois") {
    lambda <- max(meanDwell - 1, 1e-3)
    return(c(lambda = lambda, dwellSize = dwellSize))
  }
  
  if (sojournType == "shiftp") {
    lambda <- max(meanDwell - shift, 1e-3)
    return(c(lambda = lambda, shift = shift, dwellSize = dwellSize))
  }
  
  if (sojournType == "nb") {
    size <- 5
    prob <- size / (size + max(meanDwell - 1, 1e-3))
    prob <- min(max(prob, 1e-3), 1 - 1e-3)
    return(c(size = size, prob = prob, dwellSize = dwellSize))
  }
  
  if (sojournType == "shiftnb") {
    size <- 5
    targetMean <- max(meanDwell - shift, 1e-3)
    prob <- size / (size + targetMean)
    prob <- min(max(prob, 1e-3), 1 - 1e-3)
    return(c(size = size, prob = prob, shift = shift, dwellSize = dwellSize))
  }
  
  if (sojournType == "np") {
    probs <- rep(1 / dwellSize, dwellSize)
    return(probs)
  }
  
  if (sojournType == "npg") {
    # Head probabilities d(1), ..., d(m - 1) summing to less than one, plus the
    # hazard of the geometric tail. The tail hazard is chosen so that the mean of
    # the geometric part is roughly the target mean dwell time.
    headProb <- rep(0.6 / (dwellSize - 1), dwellSize - 1)
    tailProb <- min(max(1 / max(meanDwell, 1.5), 1e-3), 1 - 1e-3)
    return(c(headProb, tailProb))
  }
  
  if (sojournType == "weibull") {
    shape <- 1.5
    scale <- max(meanDwell, 1e-3)
    return(c(shape = shape, scale = scale, dwellSize = dwellSize))
  }
  
  stop("Unsupported sojournType: ", sojournType)
}


# Identical to getStateRuns() in simulate.R; kept as an internal alias since
# that is the name used throughout the initialization code below. Not
# exported (getStateRuns() is the public name for this operation).
getRunsByState <- function(states) {
  getStateRuns(states)
}


initializeMeanDwellsFromClusters <- function(x, nStates, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  clusterFit <- kmeans(
    x = matrix(x, ncol = 1),
    centers = nStates,
    nstart = 20
  )
  
  clusterOrder <- order(clusterFit$centers[, 1])
  states <- match(clusterFit$cluster, clusterOrder)
  
  stateRuns <- getRunsByState(states)
  
  meanDwells <- numeric(nStates)
  
  for (stateIndex in seq_len(nStates)) {
    dwellTimes <- stateRuns$dwellTime[stateRuns$state == stateIndex]
    
    if (length(dwellTimes) == 0L) {
      meanDwells[stateIndex] <- 2
    } else {
      meanDwells[stateIndex] <- max(mean(dwellTimes), 1.5)
    }
  }
  
  meanDwells
}


checkStateLabels <- function(stateLabels, nObs, nStates) {
  if (is.null(stateLabels)) {
    return(invisible(TRUE))
  }
  
  if (length(stateLabels) != nObs) {
    stop("stateLabels must have length equal to length(x).")
  }
  
  if (any(!stateLabels %in% seq_len(nStates))) {
    stop("stateLabels must contain only integers from 1 to nStates.")
  }
  
  invisible(TRUE)
}


initializeEmissionParametersFromLabels <- function(x, stateLabels, nStates, emissionType) {
  checkStateLabels(
    stateLabels = stateLabels,
    nObs = length(x),
    nStates = nStates
  )
  
  if (emissionType == "normal") {
    emissionParameters <- matrix(NA_real_, nrow = nStates, ncol = 2)
    
    for (stateIndex in seq_len(nStates)) {
      stateX <- x[stateLabels == stateIndex]
      
      if (length(stateX) < 2L) {
        stateX <- x
      }
      
      emissionParameters[stateIndex, 1] <- mean(stateX)
      emissionParameters[stateIndex, 2] <- max(sd(stateX), 1e-3)
    }
    
    colnames(emissionParameters) <- c("mean", "sd")
    return(emissionParameters)
  }
  
  if (emissionType == "pois") {
    emissionParameters <- matrix(NA_real_, nrow = nStates, ncol = 1)
    
    for (stateIndex in seq_len(nStates)) {
      stateX <- x[stateLabels == stateIndex]
      
      if (length(stateX) < 1L) {
        stateX <- x
      }
      
      emissionParameters[stateIndex, 1] <- max(mean(stateX), 1e-3)
    }
    
    colnames(emissionParameters) <- "lambda"
    return(emissionParameters)
  }
  
  stop("Unsupported emissionType: ", emissionType)
}


initializeEmbeddedTpmFromLabels <- function(stateLabels, nStates, pseudocount = 1) {
  tpmCounts <- matrix(pseudocount, nrow = nStates, ncol = nStates)
  diag(tpmCounts) <- 0
  
  if (length(stateLabels) >= 2L) {
    for (timeIndex in 2:length(stateLabels)) {
      fromState <- stateLabels[timeIndex - 1L]
      toState <- stateLabels[timeIndex]
      
      if (fromState != toState) {
        tpmCounts[fromState, toState] <- tpmCounts[fromState, toState] + 1
      }
    }
  }
  
  diag(tpmCounts) <- 0
  tpm <- tpmCounts / rowSums(tpmCounts)
  
  tpm
}


initializeMeanDwellsFromLabels <- function(stateLabels, nStates) {
  stateRuns <- getStateRuns(stateLabels)
  meanDwells <- numeric(nStates)
  
  for (stateIndex in seq_len(nStates)) {
    dwellTimes <- stateRuns$dwellTime[stateRuns$state == stateIndex]
    
    if (length(dwellTimes) == 0L) {
      meanDwells[stateIndex] <- 2
    } else {
      meanDwells[stateIndex] <- max(mean(dwellTimes), 1.5)
    }
  }
  
  meanDwells
}


validateEmissionParameters <- function(emissionParameters, nStates, emissionType) {
  expectedCols <- getEmissionParamCount(emissionType)
  
  if (!is.matrix(emissionParameters)) {
    stop("emissionParameters must be a matrix.")
  }
  
  if (nrow(emissionParameters) != nStates) {
    stop("emissionParameters must have nStates rows.")
  }
  
  if (ncol(emissionParameters) != expectedCols) {
    stop("emissionParameters has the wrong number of columns for emissionType.")
  }
  
  if (emissionType == "normal" && any(emissionParameters[, 2] <= 0)) {
    stop("Normal emission standard deviations must be positive.")
  }
  
  if (emissionType == "pois" && any(emissionParameters[, 1] <= 0)) {
    stop("Poisson emission rates must be positive.")
  }
  
  invisible(TRUE)
}


validateSojournParameters <- function(sojournParameters, sojournType, nStates) {
  if (!is.list(sojournParameters)) {
    stop("sojournParameters must be a list.")
  }
  
  if (length(sojournParameters) != nStates) {
    stop("sojournParameters must have one entry per state.")
  }
  
  if (length(sojournType) != nStates) {
    stop("sojournType must have one entry per state.")
  }
  
  for (stateIndex in seq_len(nStates)) {
    getDwellSize(
      sojournParameters = sojournParameters[[stateIndex]],
      sojournType = sojournType[stateIndex]
    )
  }
  
  invisible(TRUE)
}


# Perturb a starting model so that repeated starts are genuinely different.
#
# Varying only the seed is not enough: kmeans() with nstart = 20 converges to the
# same partition every time, so seed-only "multi-start" produces near-identical
# starting values and gives no protection against local optima.
jitterHsmmStart <- function(model, x, emissionScale = 0.25, dwellScale = 0.25) {
  emissionType <- model$emissionType
  emissionParameters <- model$emissionParameters
  nStates <- model$nStates
  
  if (emissionType == "normal") {
    spread <- stats::sd(x, na.rm = TRUE)
    
    emissionParameters[, 1] <- emissionParameters[, 1] +
      stats::rnorm(nStates, mean = 0, sd = emissionScale * spread)
    
    emissionParameters[, 2] <- pmax(
      emissionParameters[, 2] * exp(stats::rnorm(nStates, 0, emissionScale)),
      1e-3
    )
  }
  
  if (emissionType == "pois") {
    emissionParameters[, 1] <- pmax(
      emissionParameters[, 1] * exp(stats::rnorm(nStates, 0, emissionScale)),
      1e-3
    )
  }
  
  model$emissionParameters <- emissionParameters
  
  # Perturb the embedded transition matrix by resampling each row's off-diagonal
  # probabilities from a Dirichlet centred on the current row.
  if (nStates > 2L) {
    tpm <- model$tpm
    
    for (stateIndex in seq_len(nStates)) {
      currentRow <- tpm[stateIndex, -stateIndex]
      concentration <- pmax(currentRow, 1e-6) * 20
      gammaDraw <- stats::rgamma(length(currentRow), shape = concentration, rate = 1)
      
      if (sum(gammaDraw) > 0) {
        tpm[stateIndex, -stateIndex] <- gammaDraw / sum(gammaDraw)
      }
    }
    
    model$tpm <- tpm
  }
  
  # Perturb dwell-time parameters, leaving structural entries (shift, dwellSize) alone.
  for (stateIndex in seq_len(nStates)) {
    sojournType <- model$sojournType[stateIndex]
    parameters <- model$sojournParameters[[stateIndex]]
    
    if (sojournType %in% c("pois", "shiftp")) {
      parameters[1] <- max(parameters[1] * exp(stats::rnorm(1, 0, dwellScale)), 1e-3)
    } else if (sojournType %in% c("nb", "shiftnb", "weibull")) {
      parameters[1] <- max(parameters[1] * exp(stats::rnorm(1, 0, dwellScale)), 1e-3)
      
      if (sojournType == "weibull") {
        parameters[2] <- max(parameters[2] * exp(stats::rnorm(1, 0, dwellScale)), 1e-3)
      } else {
        logit <- stats::qlogis(min(max(parameters[2], 1e-6), 1 - 1e-6))
        parameters[2] <- stats::plogis(logit + stats::rnorm(1, 0, dwellScale))
      }
    } else if (sojournType == "geo") {
      logit <- stats::qlogis(min(max(parameters[1], 1e-6), 1 - 1e-6))
      parameters[1] <- stats::plogis(logit + stats::rnorm(1, 0, dwellScale))
    } else if (sojournType == "np") {
      weights <- stats::rgamma(length(parameters), shape = pmax(parameters, 1e-6) * 20, rate = 1)
      
      if (sum(weights) > 0) {
        parameters <- weights / sum(weights)
      }
    } else if (sojournType == "npg") {
      dwellSize <- length(parameters)
      headValues <- parameters[-dwellSize]
      headRaw <- c(1 - sum(headValues), headValues)
      weights <- stats::rgamma(length(headRaw), shape = pmax(headRaw, 1e-6) * 20, rate = 1)
      
      if (sum(weights) > 0) {
        weights <- weights / sum(weights)
        parameters[-dwellSize] <- weights[-1]
      }
      
      logit <- stats::qlogis(min(max(parameters[dwellSize], 1e-6), 1 - 1e-6))
      parameters[dwellSize] <- stats::plogis(logit + stats::rnorm(1, 0, dwellScale))
    }
    
    model$sojournParameters[[stateIndex]] <- parameters
  }
  
  model
}


#' Generate multiple perturbed starting points for HSMM fitting
#'
#' For likelihood maximization, using several starts is often preferable
#' because the likelihood can have local optima. The first start is the
#' unperturbed deterministic baseline from [initializeHsmmParameters()].
#' Every later start is a perturbation of it produced by an internal jitter
#' step: emission means are shifted by a fraction of the pooled standard
#' deviation, scale parameters are perturbed multiplicatively on the log
#' scale, and transition rows and nonparametric dwell distributions are
#' resampled from Dirichlet distributions centred on their current values.
#' Varying only the random seed is not enough for this purpose, because
#' [stats::kmeans()] with `nstart = 20` converges to the same partition
#' every time and gives no protection against local optima on its own.
#'
#' @param x Numeric vector of observations.
#' @param nStarts Integer, the number of starts to generate.
#' @param nStates Integer, the number of HSMM states.
#' @param sojournType Character vector (or single string recycled to
#'   `nStates`) naming the dwell-time family for each state.
#' @param emissionType Emission distribution, `"normal"` or `"pois"`.
#' @param dwellSize Integer vector (or single value recycled to `nStates`),
#'   the truncation/support size for each state's dwell-time distribution.
#' @param stationary Logical; whether the resulting models are intended for
#'   stationary fitting.
#' @param fitShift Logical; passed through to [buildHsmm()].
#' @param seed Optional integer seed for reproducibility.
#' @param emissionScale,dwellScale Numeric perturbation sizes controlling
#'   how far jittered starts move from the baseline.
#'
#' @return A list of `nStarts` initialization results, each with the same
#'   structure as [initializeHsmmParameters()]'s return value (`model`,
#'   `parVec`, `meanDwells`, `stateLabels`).
#'
#' @examples
#' starts <- initializeHsmmMultiStart(
#'   x = rnorm(100),
#'   nStarts = 3,
#'   nStates = 2,
#'   sojournType = c("shiftp", "np"),
#'   emissionType = "normal",
#'   dwellSize = c(10, 4),
#'   stationary = TRUE,
#'   seed = 123
#' )
#'
#' @export
initializeHsmmMultiStart <- function(
    x,
    nStarts,
    nStates,
    sojournType,
    emissionType = "normal",
    dwellSize = 10,
    stationary = TRUE,
    fitShift = FALSE,
    seed = NULL,
    emissionScale = 0.25,
    dwellScale = 0.25
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  baseStart <- initializeHsmmParameters(
    x = x,
    nStates = nStates,
    sojournType = sojournType,
    emissionType = emissionType,
    dwellSize = dwellSize,
    stationary = stationary,
    fitShift = fitShift,
    seed = NULL,
    useClusterDwells = TRUE
  )
  
  starts <- vector("list", nStarts)
  
  # The first start is the unperturbed deterministic baseline.
  starts[[1]] <- baseStart
  
  if (nStarts >= 2L) {
    for (startIndex in 2:nStarts) {
      jitteredModel <- jitterHsmmStart(
        model = baseStart$model,
        x = x,
        emissionScale = emissionScale,
        dwellScale = dwellScale
      )
      
      starts[[startIndex]] <- list(
        model = jitteredModel,
        parVec = hsmmConstrainedToWorking(
          model = jitteredModel,
          stationary = stationary
        ),
        meanDwells = baseStart$meanDwells,
        stateLabels = baseStart$stateLabels
      )
    }
  }
  
  starts
}



#' Flexible starting-value initialization for HSMM fitting
#'
#' Creates a valid constrained HSMM model and the corresponding unconstrained
#' working parameter vector used by numerical optimizers, following the
#' rule: use supplied values when available, otherwise estimate reasonable
#' starting values automatically. This lets the same function support fully
#' automatic starts, label-based starts, and manually supplied starts.
#'
#' When no additional starting values are supplied, emission parameters are
#' estimated from rough [stats::kmeans()] clustering of the observations,
#' and dwell-time parameters are initialized from either cluster-based run
#' lengths or conservative default values. Rough state labels can instead
#' be supplied through `stateLabels` (integers `1` to `nStates`, one per
#' observation) to initialize emission parameters, embedded transition
#' probabilities, and mean dwell times directly; self-transitions in
#' `stateLabels` are ignored when initializing the embedded transition
#' matrix, since it represents transitions between distinct HSMM states and
#' should have a zero diagonal. `emissionParameters`, `sojournParameters`,
#' and `tpm` can each be supplied directly to override automatic
#' initialization for that part of the model.
#'
#' @param x Numeric vector of observations.
#' @param nStates Integer, the number of HSMM states.
#' @param sojournType Character vector (or single string recycled to
#'   `nStates`) naming the dwell-time family for each state.
#' @param emissionType Emission distribution, `"normal"` or `"pois"`.
#' @param dwellSize Integer vector (or single value recycled to `nStates`),
#'   the truncation/support size for each state's dwell-time distribution.
#'   Not needed when `sojournParameters` is supplied; it is inferred from
#'   those parameters instead.
#' @param stationary Logical; whether the resulting model is intended for
#'   stationary fitting.
#' @param fitShift Logical; passed through to [buildHsmm()].
#' @param seed Optional integer seed for reproducibility of the automatic
#'   clustering steps.
#' @param useClusterDwells Logical; whether to derive mean dwell times from
#'   cluster-based run lengths (`TRUE`) or from conservative defaults based
#'   on `dwellSize` (`FALSE`) when no `stateLabels` or `sojournParameters`
#'   are supplied.
#' @param stateLabels Optional integer vector, the same length as `x`, of
#'   rough state labels in `1:nStates`.
#' @param emissionParameters Optional user-supplied emission parameter
#'   matrix, overriding automatic initialization.
#' @param sojournParameters Optional user-supplied list of dwell-time
#'   parameters, overriding automatic initialization.
#' @param tpm Optional user-supplied embedded transition matrix, overriding
#'   automatic initialization. Must have zero diagonal and rows summing to
#'   one.
#' @param delta Optional HSMM-level initial distribution, used only when
#'   `stationary = FALSE`.
#' @param pseudocount Numeric pseudocount used when estimating `tpm` from
#'   `stateLabels`.
#'
#' @return A list containing `model` (initialized constrained HSMM model),
#'   `parVec` (unconstrained working parameter vector), `meanDwells`
#'   (initial mean dwell-time estimates used to initialize dwell
#'   distributions), and `stateLabels` (the labels used for initialization,
#'   if any).
#'
#' @examples
#' init <- initializeHsmmParameters(
#'   x = rnorm(100),
#'   nStates = 2,
#'   sojournType = c("shiftp", "np"),
#'   emissionType = "normal",
#'   dwellSize = c(10, 4),
#'   stationary = TRUE,
#'   seed = 123
#' )
#' init$model
#'
#' @export
initializeHsmmParameters <- function(
    x,
    nStates,
    sojournType,
    emissionType = "normal",
    dwellSize = 10,
    stationary = TRUE,
    fitShift = FALSE,
    seed = NULL,
    useClusterDwells = TRUE,
    stateLabels = NULL,
    emissionParameters = NULL,
    sojournParameters = NULL,
    tpm = NULL,
    delta = NULL,
    pseudocount = 1
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (length(sojournType) == 1L) {
    sojournType <- rep(sojournType, nStates)
  }
  
  if (length(dwellSize) == 1L) {
    dwellSize <- rep(dwellSize, nStates)
  }
  
  if (length(sojournType) != nStates) {
    stop("sojournType must have length one or nStates.")
  }
  
  if (length(dwellSize) != nStates) {
    stop("dwellSize must have length one or nStates.")
  }
  
  checkStateLabels(
    stateLabels = stateLabels,
    nObs = length(x),
    nStates = nStates
  )
  
  if (is.null(emissionParameters)) {
    if (!is.null(stateLabels)) {
      emissionParameters <- initializeEmissionParametersFromLabels(
        x = x,
        stateLabels = stateLabels,
        nStates = nStates,
        emissionType = emissionType
      )
    } else {
      emissionParameters <- initializeEmissionParameters(
        x = x,
        nStates = nStates,
        emissionType = emissionType,
        seed = seed
      )
    }
  } else {
    validateEmissionParameters(
      emissionParameters = emissionParameters,
      nStates = nStates,
      emissionType = emissionType
    )
  }
  
  if (is.null(tpm)) {
    if (!is.null(stateLabels)) {
      tpm <- initializeEmbeddedTpmFromLabels(
        stateLabels = stateLabels,
        nStates = nStates,
        pseudocount = pseudocount
      )
    } else {
      tpm <- initializeEmbeddedTpm(nStates)
    }
  } else {
    checkTransitionMatrix(tpm, name = "tpm")
    
    if (max(abs(diag(tpm))) > 1e-8) {
      stop("tpm must be an embedded transition matrix with zero diagonal.")
    }
  }
  
  if (is.null(sojournParameters)) {
    if (!is.null(stateLabels)) {
      meanDwells <- initializeMeanDwellsFromLabels(
        stateLabels = stateLabels,
        nStates = nStates
      )
    } else if (useClusterDwells) {
      meanDwells <- initializeMeanDwellsFromClusters(
        x = x,
        nStates = nStates,
        seed = seed
      )
    } else {
      meanDwells <- pmax(2, dwellSize / 2)
    }
    
    sojournParameters <- vector("list", nStates)
    
    for (stateIndex in seq_len(nStates)) {
      sojournParameters[[stateIndex]] <- initializeSojournOne(
        sojournType = sojournType[stateIndex],
        dwellSize = dwellSize[stateIndex],
        meanDwell = meanDwells[stateIndex],
        shift = 1,
        fitShift = fitShift
      )
    }
  } else {
    validateSojournParameters(
      sojournParameters = sojournParameters,
      sojournType = sojournType,
      nStates = nStates
    )
    
    meanDwells <- rep(NA_real_, nStates)
  }
  
  if (stationary) {
    deltaForModel <- NULL
  } else {
    if (is.null(delta)) {
      if (!is.null(stateLabels)) {
        delta <- tabulate(stateLabels[1], nbins = nStates)
        delta <- pmax(delta, 1e-3)
        delta <- delta / sum(delta)
      } else {
        delta <- rep(1 / nStates, nStates)
      }
    }
    
    deltaForModel <- delta
  }
  
  model <- buildHsmm(
    tpm = tpm,
    sojournParameters = sojournParameters,
    sojournType = sojournType,
    emissionParameters = emissionParameters,
    emissionType = emissionType,
    delta = deltaForModel,
    stationary = stationary,
    fitShift = fitShift
  )
  
  parVec <- hsmmConstrainedToWorking(
    model = model,
    stationary = stationary
  )
  
  list(
    model = model,
    parVec = parVec,
    meanDwells = meanDwells,
    stateLabels = stateLabels
  )
}
