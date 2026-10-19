# Emission-family registry for the EM fitter.
#
# Every emission distribution is a small object with the same interface, so
# em.R never branches on emissionType:
#
#   prepare(x, control)                     -> data reduced to what the density needs
#   logDensity(prepared, params, control)   -> nObs x nStates log-density matrix
#   update(prepared, weights, params, control, fitControl) -> new params (M-step)
#   simulate(states, params, control)       -> data for an HSMM-level state path
#   initialize(prepared, nStates, control, method) -> list(params, labels)
#   exactMStep(control)                     -> TRUE if update() maximizes Q exactly
#
# `params` is always an nStates x p matrix (one row per HSMM state), which is
# what buildHsmm() and hsmmToHmm() expect.
#
# Scalar families: "normal", "pois" (x is a numeric vector).
#
# Functional (topological HMM) families, Kashlak, Loliencar and Heo (2023).
# x is an nObs x L matrix (one curve per row) or an nObs x L x d array
# (d-dimensional curves, coordinates independent), observed on
# control$grid (default seq(0, 1, length.out = L)):
#
#   "bmwd"   Brownian motion with drift, optionally fractional (control$hurst).
#            Onsager-Machlup functional, constant drift c (paper Sec. 4.2-4.3):
#              log b_j(O) = c_j S(O) - c_j^2 Q / 2,
#              S(O) = int g dO,  Q = int g^2 dtau,  g(tau) = kappa tau^ups,
#              ups = 1/2 - hurst,  kappa = Gamma(ups + 1) / Gamma(2 ups + 1),
#            evaluated exactly at the piecewise-linear interpolant of O.
#            hurst = 1/2 gives S = O(1) - O(0), Q = 1 (Girsanov). The M-step
#            c = sum w S / (Q sum w) is the paper's closed form, Gamma factors
#            included. "BM2d" in thmm1.R is the d = 2 case.
#   "ou"     Ornstein-Uhlenbeck dO = (b0 - b1 O) dtau + dW in the paper's
#            (b0, b1) = (c mu, c) parametrization (Sec. 5.2). The functional
#              log b = b0 Delta - b1 I - (b0^2 span - 2 b0 b1 M1 + b1^2 M2) / 2
#            (Delta = sum dO, I = sum O dO (Ito, left point), M1 = int O,
#            M2 = int O^2) is the paper's -1/2 int (O' - b)^2 + b1/2 written
#            with an Ito integral; it is a concave quadratic in (b0, b1), so the
#            M-step is an exact 2 x 2 linear solve (with b1 >= 0), not optim().
#   "nonpar" Nonparametric state mean curves, log b_j = -scale/2 ||O - h_j||^2
#            with norm "L2", "W21", "W22" (sums of squared values / first /
#            second differences on the grid, exactly as thmm1.R) or "sup"
#            (max |O - h_j|, unsquared, as thmm1.R). The M-step is the weighted
#            mean curve: exact for the three quadratic norms, a heuristic for
#            "sup" (exactMStep = FALSE).
#
# Where thmm1.R departs from the paper's functionals (BMWD uses
# (sum |dO - c/n|)^2; OU uses the square of the summed residual; fractional BM
# uses the BM emission in the E-step but the fBM drift in the M-step) these
# families follow the paper, so that each M-step maximizes the same functional
# its E-step evaluates.


`%||%` <- function(a, b) if (is.null(a)) b else a


getEmissionFamily <- function(emissionType) {
  switch(
    emissionType,
    normal = normalEmissionFamily(),
    pois = poissonEmissionFamily(),
    bmwd = bmwdEmissionFamily(),
    ou = ouEmissionFamily(),
    nonpar = nonparEmissionFamily(),
    stop("Unsupported emissionType: ", emissionType)
  )
}


isFunctionalEmission <- function(emissionType) {
  getEmissionFamily(emissionType)$isFunctional
}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Hard (one-hot) weights from integer labels: lets update() double as the
# "fit each cluster separately" step of initialization.
labelsToWeights <- function(labels, nStates) {
  weights <- matrix(0, nrow = length(labels), ncol = nStates)
  observed <- !is.na(labels)
  weights[cbind(which(observed), labels[observed])] <- 1
  weights
}


# Shared initialization: cluster per-observation features, fit each cluster
# with the family's own M-step, or (method = "random") fit single random
# observations, as thmm1.R's autoParam() does.
initializeFromFeatures <- function(family, prepared, features, nStates, control, method, orderBy) {
  observed <- which(!prepared$missing)
  labels <- rep(NA_integer_, length(prepared$missing))

  if (method == "kmeans") {
    features <- as.matrix(features)[observed, , drop = FALSE]
    clusterFit <- kmeans(features, centers = nStates, nstart = 20)
    labels[observed] <- clusterFit$cluster
  } else if (method == "random") {
    chosen <- sample(observed, nStates)
    labels[chosen] <- seq_len(nStates)
  } else {
    stop("method must be \"kmeans\" or \"random\".")
  }

  params <- family$update(
    prepared,
    labelsToWeights(labels, nStates),
    NULL,
    control,
    list(minCount = 0, minSd = 1e-3)
  )

  stateOrder <- order(orderBy(params))
  params <- params[stateOrder, , drop = FALSE]

  # Labels for the transition / dwell starting values: most likely state of each
  # observation under the initial emission parameters, ignoring dynamics.
  logDensity <- family$logDensity(prepared, params, control)
  labels <- max.col(logDensity, ties.method = "first")

  list(params = params, labels = labels)
}


weightTotals <- function(weights, missing) {
  weights[missing, ] <- 0
  list(weights = weights, totals = colSums(weights))
}


# ---------------------------------------------------------------------------
# Scalar families
# ---------------------------------------------------------------------------

normalEmissionFamily <- function() {
  list(
    name = "normal",
    isFunctional = FALSE,
    exactMStep = function(control) TRUE,

    prepare = function(x, control) {
      if (!is.null(dim(x))) stop("Normal emissions need a numeric vector.")
      list(x = as.numeric(x), missing = is.na(x), nObs = length(x))
    },

    logDensity = function(prepared, params, control) {
      if (any(params[, 2] <= 0)) stop("Normal emission standard deviations must be positive.")
      out <- matrix(0, nrow = prepared$nObs, ncol = nrow(params))
      obs <- !prepared$missing
      for (j in seq_len(nrow(params))) {
        out[obs, j] <- dnorm(prepared$x[obs], params[j, 1], params[j, 2], log = TRUE)
      }
      out
    },

    update = function(prepared, weights, params, control, fitControl) {
      wt <- weightTotals(weights, prepared$missing)
      x <- ifelse(prepared$missing, 0, prepared$x)
      nStates <- ncol(weights)
      updated <- params %||% matrix(c(0, 1), nrow = nStates, ncol = 2, byrow = TRUE)
      for (j in seq_len(nStates)) {
        if (wt$totals[j] <= fitControl$minCount) next
        w <- wt$weights[, j]
        mean <- sum(w * x) / wt$totals[j]
        variance <- sum(w * (x - mean)^2) / wt$totals[j]
        updated[j, ] <- c(mean, max(sqrt(variance), fitControl$minSd))
      }
      colnames(updated) <- c("mean", "sd")
      updated
    },

    simulate = NULL,   # simulateHsmm() keeps its original scalar path

    initialize = function(prepared, nStates, control, method = "kmeans") {
      family <- normalEmissionFamily()
      initializeFromFeatures(
        family, prepared, matrix(prepared$x, ncol = 1), nStates, control, method,
        orderBy = function(p) p[, 1]
      )
    }
  )
}


poissonEmissionFamily <- function() {
  list(
    name = "pois",
    isFunctional = FALSE,
    exactMStep = function(control) TRUE,

    prepare = function(x, control) {
      if (!is.null(dim(x))) stop("Poisson emissions need a numeric vector.")
      list(x = as.numeric(x), missing = is.na(x), nObs = length(x))
    },

    logDensity = function(prepared, params, control) {
      out <- matrix(0, nrow = prepared$nObs, ncol = nrow(params))
      obs <- !prepared$missing
      for (j in seq_len(nrow(params))) {
        out[obs, j] <- dpois(prepared$x[obs], params[j, 1], log = TRUE)
      }
      out
    },

    update = function(prepared, weights, params, control, fitControl) {
      wt <- weightTotals(weights, prepared$missing)
      x <- ifelse(prepared$missing, 0, prepared$x)
      nStates <- ncol(weights)
      updated <- params %||% matrix(1, nrow = nStates, ncol = 1)
      for (j in seq_len(nStates)) {
        if (wt$totals[j] <= fitControl$minCount) next
        updated[j, 1] <- max(sum(wt$weights[, j] * x) / wt$totals[j], 1e-8)
      }
      colnames(updated) <- "lambda"
      updated
    },

    simulate = NULL,

    initialize = function(prepared, nStates, control, method = "kmeans") {
      family <- poissonEmissionFamily()
      initializeFromFeatures(
        family, prepared, matrix(prepared$x, ncol = 1), nStates, control, method,
        orderBy = function(p) p[, 1]
      )
    }
  )
}


# ---------------------------------------------------------------------------
# Curve helpers
# ---------------------------------------------------------------------------

asCurveArray <- function(x) {
  if (is.null(dim(x))) {
    stop("Functional emissions need curves: an nObs x L matrix or an nObs x L x d array.")
  }

  if (length(dim(x)) == 2L) {
    x <- array(x, dim = c(dim(x), 1L))
  }

  if (length(dim(x)) != 3L) {
    stop("Curve data must be a matrix (nObs x L) or a 3-d array (nObs x L x d).")
  }

  if (dim(x)[2] < 3L) {
    stop("Each curve needs at least 3 grid points.")
  }

  x
}


getCurveGrid <- function(control, nPoints) {
  grid <- control$grid %||% seq(0, 1, length.out = nPoints)

  if (length(grid) != nPoints) {
    stop("emissionControl$grid has length ", length(grid), " but curves have ", nPoints, " points.")
  }

  if (any(diff(grid) <= 0)) {
    stop("emissionControl$grid must be strictly increasing.")
  }

  grid
}


# A curve is missing only if every value is NA; partially observed curves
# would need a modified functional and are rejected.
findMissingCurves <- function(curves) {
  perCurve <- prod(dim(curves)[-1])
  naCount <- rowSums(matrix(is.na(curves), nrow = dim(curves)[1]))

  if (any(naCount > 0 & naCount < perCurve)) {
    stop("Some curves are partially missing; only entirely missing curves are supported.")
  }

  naCount == perCurve
}


# Return simulated curves as a matrix when d = 1, an array otherwise.
dropSingleCoordinate <- function(curves) {
  if (dim(curves)[3] == 1L) {
    matrix(curves[, , 1], nrow = dim(curves)[1])
  } else {
    curves
  }
}


getSimulationGrid <- function(control) {
  control$grid %||% seq(0, 1, length.out = control$nPoints %||% 100)
}


# ---------------------------------------------------------------------------
# Brownian motion with drift (standard or fractional)
# ---------------------------------------------------------------------------

bmwdDriftKernel <- function(grid, hurst) {
  if (!is.numeric(hurst) || hurst <= 0 || hurst >= 1) {
    stop("emissionControl$hurst must lie in (0, 1).")
  }

  # The functional is evaluated exactly at the piecewise-linear interpolant of
  # the observed curve (the paper asks for observations to be projected into
  # the Cameron-Martin space before differentiating): with elapsed time
  # s = tau - tau_1,
  #   int g dO  = sum_k gbar_k (O_{k+1} - O_k),  gbar_k = cell average of g,
  #   int g^2   = kappa^2 s_L^(2 ups + 1) / (2 ups + 1).
  # Cell averages and the exact integral matter for hurst > 1/2, where
  # g ~ tau^ups is singular at 0; for hurst = 1/2 both reduce to g = 1, Q = span.
  upsilon <- 0.5 - hurst
  kappa <- gamma(upsilon + 1) / gamma(2 * upsilon + 1)
  elapsed <- grid - grid[1]
  width <- diff(elapsed)
  power <- upsilon + 1
  kernel <- kappa * (elapsed[-1]^power - elapsed[-length(elapsed)]^power) / (power * width)
  Q <- kappa^2 * elapsed[length(elapsed)]^(2 * upsilon + 1) / (2 * upsilon + 1)

  list(kernel = kernel, Q = Q)
}


bmwdEmissionFamily <- function() {
  family <- list(
    name = "bmwd",
    isFunctional = TRUE,
    exactMStep = function(control) TRUE
  )

  family$prepare <- function(x, control) {
    curves <- asCurveArray(x)
    missing <- findMissingCurves(curves)
    nPoints <- dim(curves)[2]
    grid <- getCurveGrid(control, nPoints)
    drift <- bmwdDriftKernel(grid, control$hurst %||% 0.5)

    # S[t, d] = sum_k g_k (O_{k+1} - O_k): the only data summary the
    # functional depends on, so the E- and M-steps never touch the curves.
    S <- matrix(0, nrow = dim(curves)[1], ncol = dim(curves)[3])
    for (d in seq_len(dim(curves)[3])) {
      increments <- curves[, -1, d, drop = FALSE] - curves[, -nPoints, d, drop = FALSE]
      dim(increments) <- c(dim(curves)[1], nPoints - 1L)
      S[, d] <- increments %*% drift$kernel
    }
    S[missing, ] <- 0

    list(S = S, Q = drift$Q, missing = missing, nObs = dim(curves)[1])
  }

  family$logDensity <- function(prepared, params, control) {
    out <- prepared$S %*% t(params) -
      matrix(0.5 * prepared$Q * rowSums(params^2), nrow = prepared$nObs, ncol = nrow(params), byrow = TRUE)
    out[prepared$missing, ] <- 0
    out
  }

  family$update <- function(prepared, weights, params, control, fitControl) {
    wt <- weightTotals(weights, prepared$missing)
    nStates <- ncol(weights)
    d <- ncol(prepared$S)
    updated <- params %||% matrix(0, nrow = nStates, ncol = d)
    estimate <- crossprod(wt$weights, prepared$S) / (prepared$Q * wt$totals)
    keep <- wt$totals > fitControl$minCount
    updated[keep, ] <- estimate[keep, , drop = FALSE]
    colnames(updated) <- if (d == 1L) "drift" else paste0("drift", seq_len(d))
    updated
  }

  family$simulate <- function(states, params, control) {
    grid <- getSimulationGrid(control)
    hurst <- control$hurst %||% 0.5
    nObs <- length(states)
    nPoints <- length(grid)
    d <- ncol(params)
    elapsed <- grid - grid[1]
    curves <- array(control$x0 %||% 0, dim = c(nObs, nPoints, d))

    if (hurst == 0.5) {
      width <- diff(grid)
      for (k in seq_len(d)) {
        noise <- matrix(rnorm(nObs * (nPoints - 1L)), nObs) * rep(sqrt(width), each = nObs)
        increments <- noise + outer(params[states, k], width)
        curves[, , k] <- curves[, , k] + cbind(0, t(apply(increments, 1, cumsum)))
      }
    } else {
      # fBM on the grid plus a linear drift c * tau, as in the paper (Sec. 5.3).
      s <- elapsed[-1]
      covariance <- 0.5 * (outer(s^(2 * hurst), s^(2 * hurst), "+") - abs(outer(s, s, "-"))^(2 * hurst))
      eig <- eigen(covariance, symmetric = TRUE)
      root <- eig$vectors %*% (sqrt(pmax(eig$values, 0)) * t(eig$vectors))
      for (k in seq_len(d)) {
        noise <- matrix(rnorm(nObs * (nPoints - 1L)), nObs) %*% root
        curves[, -1, k] <- curves[, -1, k] + noise
        curves[, , k] <- curves[, , k] + outer(params[states, k], elapsed)
      }
    }

    dropSingleCoordinate(curves)
  }

  family$initialize <- function(prepared, nStates, control, method = "kmeans") {
    initializeFromFeatures(
      family, prepared, prepared$S / prepared$Q, nStates, control, method,
      orderBy = function(p) p[, 1]
    )
  }

  family
}


# ---------------------------------------------------------------------------
# Ornstein-Uhlenbeck
# ---------------------------------------------------------------------------

# Exact maximizer of the concave quadratic
#   b0 SD - b1 SI - (b0^2 Sspan - 2 b0 b1 SM1 + b1^2 SM2) / 2   subject to b1 >= 0.
solveOuQuadratic <- function(SD, SI, SM1, SM2, Sspan) {
  determinant <- Sspan * SM2 - SM1^2
  rate <- (SM1 * SD - Sspan * SI) / determinant
  intercept <- (SD * SM2 - SM1 * SI) / determinant

  # Boundary of b1 >= 0: with b1 = 0 the functional is Brownian motion with drift.
  atBoundary <- !is.finite(rate) | rate < 0
  rate[atBoundary] <- 0
  intercept[atBoundary] <- SD[atBoundary] / Sspan[atBoundary]

  list(intercept = intercept, rate = rate)
}


ouEmissionFamily <- function() {
  family <- list(
    name = "ou",
    isFunctional = TRUE,
    exactMStep = function(control) TRUE
  )

  family$prepare <- function(x, control) {
    curves <- asCurveArray(x)
    missing <- findMissingCurves(curves)
    nObs <- dim(curves)[1]
    nPoints <- dim(curves)[2]
    d <- dim(curves)[3]
    width <- diff(getCurveGrid(control, nPoints))

    stats <- lapply(c("Delta", "I", "M1", "M2"), function(n) matrix(0, nObs, d))
    names(stats) <- c("Delta", "I", "M1", "M2")

    for (k in seq_len(d)) {
      path <- matrix(curves[, , k], nrow = nObs)
      left <- path[, -nPoints, drop = FALSE]
      increments <- path[, -1, drop = FALSE] - left
      stats$Delta[, k] <- rowSums(increments)
      stats$I[, k] <- rowSums(left * increments)      # Ito (left-point) integral
      stats$M1[, k] <- left %*% width
      stats$M2[, k] <- left^2 %*% width
    }

    for (n in names(stats)) stats[[n]][missing, ] <- 0

    c(stats, list(span = sum(width), d = d, missing = missing, nObs = nObs))
  }

  family$logDensity <- function(prepared, params, control) {
    d <- prepared$d
    out <- matrix(0, nrow = prepared$nObs, ncol = nrow(params))

    for (k in seq_len(d)) {
      b0 <- params[, k]
      b1 <- params[, d + k]
      out <- out + outer(prepared$Delta[, k], b0) - outer(prepared$I[, k], b1) -
        0.5 * (matrix(prepared$span * b0^2, prepared$nObs, length(b0), byrow = TRUE) -
                 2 * outer(prepared$M1[, k], b0 * b1) + outer(prepared$M2[, k], b1^2))
    }

    out[prepared$missing, ] <- 0
    out
  }

  family$update <- function(prepared, weights, params, control, fitControl) {
    wt <- weightTotals(weights, prepared$missing)
    nStates <- ncol(weights)
    d <- prepared$d
    updated <- params %||% matrix(0, nrow = nStates, ncol = 2 * d)
    keep <- wt$totals > fitControl$minCount

    for (k in seq_len(d)) {
      solution <- solveOuQuadratic(
        SD = as.vector(crossprod(wt$weights, prepared$Delta[, k])),
        SI = as.vector(crossprod(wt$weights, prepared$I[, k])),
        SM1 = as.vector(crossprod(wt$weights, prepared$M1[, k])),
        SM2 = as.vector(crossprod(wt$weights, prepared$M2[, k])),
        Sspan = prepared$span * wt$totals
      )
      updated[keep, k] <- solution$intercept[keep]
      updated[keep, d + k] <- solution$rate[keep]
    }

    colnames(updated) <- if (d == 1L) c("intercept", "rate") else
      c(paste0("intercept", seq_len(d)), paste0("rate", seq_len(d)))
    updated
  }

  family$simulate <- function(states, params, control) {
    grid <- getSimulationGrid(control)
    width <- diff(grid)
    nObs <- length(states)
    nPoints <- length(grid)
    d <- ncol(params) / 2
    curves <- array(control$x0 %||% 0, dim = c(nObs, nPoints, d))

    for (k in seq_len(d)) {
      b0 <- params[states, k]
      b1 <- params[states, d + k]
      for (step in seq_len(nPoints - 1L)) {          # Euler scheme, as thmm1.R
        current <- curves[, step, k]
        curves[, step + 1L, k] <- current + (b0 - b1 * current) * width[step] +
          sqrt(width[step]) * rnorm(nObs)
      }
    }

    dropSingleCoordinate(curves)
  }

  family$initialize <- function(prepared, nStates, control, method = "kmeans") {
    # Per-curve fits are noisy in the rate, so cluster on standardized
    # (intercept, rate) and then refit each cluster with pooled statistics.
    d <- prepared$d
    perCurve <- lapply(seq_len(d), function(k) solveOuQuadratic(
      prepared$Delta[, k], prepared$I[, k], prepared$M1[, k], prepared$M2[, k],
      rep(prepared$span, prepared$nObs)
    ))
    features <- do.call(cbind, lapply(perCurve, function(s) {
      rate <- pmin(s$rate, quantile(s$rate, 0.95, na.rm = TRUE))
      cbind(s$intercept, rate)
    }))
    features <- scale(features)
    features[!is.finite(features)] <- 0

    initializeFromFeatures(
      family, prepared, features, nStates, control, method,
      orderBy = function(p) p[, 1] / pmax(p[, d + 1], 1e-8)   # order by mean level
    )
  }

  family
}


ouMeanLevel <- function(params) {
  d <- ncol(params) / 2
  params[, seq_len(d), drop = FALSE] / params[, d + seq_len(d), drop = FALSE]
}


# ---------------------------------------------------------------------------
# Nonparametric mean curves
# ---------------------------------------------------------------------------

# Difference operator for one coordinate: identity (L2), first differences
# (W21) or second differences (W22).
makeNormOperator <- function(norm, nPoints) {
  switch(
    norm,
    L2 = diag(nPoints),
    W21 = diff(diag(nPoints)),
    W22 = diff(diag(nPoints), differences = 2),
    stop("Unsupported norm: ", norm)
  )
}


nonparEmissionFamily <- function() {
  family <- list(
    name = "nonpar",
    isFunctional = TRUE,
    exactMStep = function(control) !identical(control$norm %||% "L2", "sup")
  )

  family$prepare <- function(x, control) {
    curves <- asCurveArray(x)
    missing <- findMissingCurves(curves)
    nObs <- dim(curves)[1]
    nPoints <- dim(curves)[2]
    d <- dim(curves)[3]
    norm <- control$norm %||% "L2"

    # Rows are c(curve[, 1], ..., curve[, d]), matching the stored mean curves.
    flat <- matrix(curves, nrow = nObs)
    flat[missing, ] <- 0

    prepared <- list(flat = flat, norm = norm, nPoints = nPoints, d = d,
                     scale = control$scale %||% 1, missing = missing, nObs = nObs)

    if (norm != "sup") {
      # ||P(O - h)||^2 = ||PO||^2 - 2 <PO, Ph> + ||Ph||^2 with P block-diagonal
      # over coordinates; PO is computed once per fit.
      operator <- makeNormOperator(norm, nPoints)
      prepared$operator <- kronecker(diag(d), operator)
      prepared$projected <- flat %*% t(prepared$operator)
      prepared$projectedSq <- rowSums(prepared$projected^2)
    }

    prepared
  }

  family$logDensity <- function(prepared, params, control) {
    nStates <- nrow(params)

    if (prepared$norm == "sup") {
      distance <- matrix(0, prepared$nObs, nStates)
      for (j in seq_len(nStates)) {
        absDiff <- abs(sweep(prepared$flat, 2, params[j, ]))
        for (k in seq_len(prepared$d)) {
          cols <- (k - 1L) * prepared$nPoints + seq_len(prepared$nPoints)
          distance[, j] <- distance[, j] + apply(absDiff[, cols, drop = FALSE], 1, max)
        }
      }
    } else {
      projectedMeans <- params %*% t(prepared$operator)
      distance <- prepared$projectedSq - 2 * prepared$projected %*% t(projectedMeans) +
        matrix(rowSums(projectedMeans^2), prepared$nObs, nStates, byrow = TRUE)
      distance <- pmax(distance, 0)
    }

    out <- -0.5 * prepared$scale * distance
    out[prepared$missing, ] <- 0
    out
  }

  family$update <- function(prepared, weights, params, control, fitControl) {
    wt <- weightTotals(weights, prepared$missing)
    updated <- params %||% matrix(0, ncol(weights), ncol(prepared$flat))
    keep <- wt$totals > fitControl$minCount
    means <- crossprod(wt$weights, prepared$flat) / wt$totals
    updated[keep, ] <- means[keep, , drop = FALSE]
    updated
  }

  family$simulate <- function(states, params, control) {
    # Mean curve plus the truncated Karhunen-Loeve expansion of a Brownian bridge
    # used in the paper (Sec. 5.4) and in thmm1.R's rSmoothError().
    grid <- getSimulationGrid(control)
    nPoints <- length(grid)
    d <- ncol(params) / nPoints
    nTerms <- control$noiseTerms %||% 16
    noiseSd <- control$noiseSd %||% 0.4
    basis <- sqrt(2) * sin(outer(seq_len(nTerms), pi * grid)) / (seq_len(nTerms) * pi)
    nObs <- length(states)
    flat <- params[states, , drop = FALSE]

    for (k in seq_len(d)) {
      cols <- (k - 1L) * nPoints + seq_len(nPoints)
      coefficients <- matrix(rnorm(nObs * nTerms, sd = noiseSd), nObs, nTerms)
      flat[, cols] <- flat[, cols] + coefficients %*% basis
    }

    dropSingleCoordinate(array(flat, dim = c(nObs, nPoints, d)))
  }

  family$initialize <- function(prepared, nStates, control, method = "kmeans") {
    # Functional k-means on the grid (L2), as recommended in the paper (Sec. 5.4).
    initializeFromFeatures(
      family, prepared, prepared$flat, nStates, control, method,
      orderBy = function(p) p[, 1]
    )
  }

  family
}


# ---------------------------------------------------------------------------
# Starting values for fitHsmmEM()
# ---------------------------------------------------------------------------

#' Starting values for EM fitting, for any emission family
#'
#' Builds a starting HSMM for [fitHsmmEM()] from the data: emission
#' parameters from the family's own initializer (k-means on per-observation
#' summaries followed by a per-cluster fit, or single random observations as
#' in thmm1.R's `autoParam()`), then the embedded transition matrix and mean
#' dwell times from the most likely state of each observation.
#'
#' @param x Data: a numeric vector (scalar families) or an nObs x L matrix /
#'   nObs x L x d array of curves (functional families).
#' @param nStates Number of HSMM states.
#' @param sojournType Dwell-time family per state (recycled).
#' @param emissionType One of `"normal"`, `"pois"`, `"bmwd"`, `"ou"`,
#'   `"nonpar"`.
#' @param emissionControl Family settings, e.g. `list(grid, hurst)` for
#'   `"bmwd"` or `list(norm, scale)` for `"nonpar"`.
#' @param dwellSize Aggregate size per state (recycled).
#' @param stationary Logical, as in [fitHsmmEM()].
#' @param method `"kmeans"` or `"random"`.
#' @param seed Optional seed.
#'
#' @return A list with `model` (the starting HSMM), `stateLabels`, and
#'   `meanDwells`.
#' @export
initializeHsmmEM <- function(
    x,
    nStates,
    sojournType,
    emissionType,
    emissionControl = NULL,
    dwellSize = 10,
    stationary = FALSE,
    method = c("kmeans", "random"),
    seed = NULL
) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  sojournType <- rep_len(sojournType, nStates)
  dwellSize <- rep_len(dwellSize, nStates)
  family <- getEmissionFamily(emissionType)
  prepared <- family$prepare(x, emissionControl)
  start <- family$initialize(prepared, nStates, emissionControl, method)

  labels <- start$labels
  tpm <- initializeEmbeddedTpmFromLabels(labels, nStates)
  meanDwells <- initializeMeanDwellsFromLabels(labels, nStates)

  sojournParameters <- lapply(seq_len(nStates), function(i) {
    initializeSojournOne(sojournType[i], dwellSize = dwellSize[i], meanDwell = meanDwells[i])
  })

  delta <- tabulate(labels[1], nbins = nStates) + 0.1
  delta <- delta / sum(delta)

  model <- buildHsmm(
    tpm = tpm,
    sojournParameters = sojournParameters,
    sojournType = sojournType,
    emissionParameters = start$params,
    emissionType = emissionType,
    emissionControl = emissionControl,
    delta = delta,
    stationary = stationary
  )

  list(model = model, stateLabels = labels, meanDwells = meanDwells)
}
