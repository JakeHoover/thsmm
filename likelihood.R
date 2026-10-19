#' Compute the state-dependent emission probability matrix
#'
#' Evaluates the emission density (or pmf) of every observation under every
#' expanded HMM state.
#'
#' @param x Numeric vector of observations. `NA` entries are treated as
#'   missing: they are given a state-dependent likelihood of one in every
#'   state, which makes the forward recursion skip them, the standard HMM
#'   treatment of missing data — series with gaps can be passed in directly.
#' @param stateParameters Numeric matrix of expanded state-dependent
#'   emission parameters, one row per expanded HMM state (e.g.
#'   `hmmModel$stateParameters` from [hsmmToHmm()]).
#' @param emissionType One of `"normal"` (columns `mean`, `sd`; the second
#'   parameter is the standard deviation, not the variance) or `"pois"`
#'   (column `lambda`).
#' @param minDensity Numeric floor applied to every entry, to keep the
#'   likelihood strictly positive.
#'
#' @return A numeric matrix with `length(x)` rows and `nrow(stateParameters)`
#'   columns.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' hmmModel <- hsmmToHmm(testModel)
#' emissionMatrix <- makeEmissionMatrix(
#'   x = c(0.1, 4.2, NA, -0.3),
#'   stateParameters = hmmModel$stateParameters,
#'   emissionType = hmmModel$emissionType
#' )
#'
#' @export
makeEmissionMatrix <- function(x, stateParameters, emissionType, minDensity = 1e-300) {
  nObs <- length(x)
  nStates <- nrow(stateParameters)
  
  if (nObs < 1L) {
    stop("x must contain at least one observation.")
  }
  
  emissionMatrix <- matrix(NA_real_, nrow = nObs, ncol = nStates)
  
  missingIndex <- is.na(x)
  
  if (all(missingIndex)) {
    stop("x contains no observed values.")
  }
  
  if (emissionType == "normal") {
    mu <- stateParameters[, 1]
    sd <- stateParameters[, 2]
    
    if (any(sd <= 0) || any(!is.finite(sd))) {
      stop("Normal emission standard deviations must be positive and finite.")
    }
    
    for (stateIndex in seq_len(nStates)) {
      emissionMatrix[, stateIndex] <- dnorm(
        x = x,
        mean = mu[stateIndex],
        sd = sd[stateIndex]
      )
    }
  } else if (emissionType == "pois") {
    lambda <- stateParameters[, 1]
    
    if (any(lambda <= 0) || any(!is.finite(lambda))) {
      stop("Poisson emission rates must be positive and finite.")
    }
    
    for (stateIndex in seq_len(nStates)) {
      emissionMatrix[, stateIndex] <- dpois(
        x = x,
        lambda = lambda[stateIndex]
      )
    }
  } else {
    stop("Unsupported emissionType: ", emissionType)
  }
  
  # Missing observations contribute no information: setting their state-dependent
  # likelihood to one makes the forward recursion skip them, which is the standard
  # HMM treatment of missing data.
  if (any(missingIndex)) {
    emissionMatrix[missingIndex, ] <- 1
  }
  
  if (any(!is.finite(emissionMatrix))) {
    stop("Emission matrix contains non-finite values.")
  }
  
  pmax(emissionMatrix, minDensity)
}


#' Scaled forward recursion for an HMM
#'
#' Runs the numerically stable scaled-forward algorithm for the expanded
#' HMM, the single implementation the whole package uses to evaluate the
#' likelihood. This function sits on the hot path of every likelihood
#' evaluation during fitting; see [hsmmNegativeLogLik()] for the full
#' pipeline that builds its inputs from a working parameter vector.
#'
#' @param delta Numeric vector, the expanded HMM initial distribution.
#' @param gamma Numeric matrix, the expanded HMM transition matrix.
#' @param emissionMatrix Numeric matrix of state-dependent likelihoods, as
#'   returned by [makeEmissionMatrix()].
#' @param returnAlpha Logical; whether to include the scaled forward
#'   probabilities in the result.
#' @param tol Numeric tolerance used when clamping and renormalizing `delta`
#'   and `gamma`.
#'
#' @return A list containing `logLik`, `negativeLogLik`, `scales` (the
#'   forward scaling constants), and, if `returnAlpha = TRUE`, `alpha` (the
#'   scaled forward probabilities; rows should sum to one). If the recursion
#'   collapses (a scale factor of zero or a non-finite value), the same
#'   fields are returned with `logLik = -Inf`, `negativeLogLik = Inf`,
#'   `scales` all `NA`, and `alpha = NULL`, so the returned shape is always
#'   the same.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' hmmModel <- hsmmToHmm(testModel)
#' emissionMatrix <- makeEmissionMatrix(
#'   x = rnorm(20),
#'   stateParameters = hmmModel$stateParameters,
#'   emissionType = hmmModel$emissionType
#' )
#' forwardResult <- scaledForward(
#'   delta = hmmModel$delta,
#'   gamma = hmmModel$gamma,
#'   emissionMatrix = emissionMatrix,
#'   returnAlpha = TRUE
#' )
#' max(abs(rowSums(forwardResult$alpha) - 1))
#'
#' @export
scaledForward <- function(delta, gamma, emissionMatrix, returnAlpha = FALSE, tol = 1e-12) {
  nObs <- nrow(emissionMatrix)
  nStates <- ncol(emissionMatrix)
  
  if (length(delta) != nStates) {
    stop("delta length does not match number of HMM states.")
  }
  
  if (nrow(gamma) != nStates || ncol(gamma) != nStates) {
    stop("gamma dimensions do not match number of HMM states.")
  }
  
  if (any(delta < -tol)) {
    stop("delta contains negative values.")
  }
  
  if (any(gamma < -tol)) {
    stop("gamma contains negative values.")
  }
  
  delta <- pmax(delta, 0)
  delta <- delta / sum(delta)
  
  gamma <- pmax(gamma, 0)
  gamma <- gamma / rowSums(gamma)
  
  alpha <- matrix(NA_real_, nrow = nObs, ncol = nStates)
  scales <- numeric(nObs)
  
  failedResult <- list(
    logLik = -Inf,
    negativeLogLik = Inf,
    scales = rep(NA_real_, nObs),
    alpha = NULL
  )
  
  alphaRaw <- delta * emissionMatrix[1, ]
  scales[1] <- sum(alphaRaw)
  
  if (!is.finite(scales[1]) || scales[1] <= 0) {
    return(failedResult)
  }
  
  alpha[1, ] <- alphaRaw / scales[1]
  
  if (nObs >= 2L) {
    for (timeIndex in 2:nObs) {
      alphaRaw <- as.numeric(alpha[timeIndex - 1, ] %*% gamma) * emissionMatrix[timeIndex, ]
      scales[timeIndex] <- sum(alphaRaw)
      
      if (!is.finite(scales[timeIndex]) || scales[timeIndex] <= 0) {
        return(failedResult)
      }
      
      alpha[timeIndex, ] <- alphaRaw / scales[timeIndex]
    }
  }
  
  logLik <- sum(log(scales))
  
  output <- list(
    logLik = logLik,
    negativeLogLik = -logLik,
    scales = scales
  )
  
  if (returnAlpha) {
    output$alpha <- alpha
  }
  
  output
}


#' Negative log-likelihood of an HSMM, for optimization
#'
#' The full objective function for fitting an HSMM by maximum likelihood.
#' Runs the pipeline `parVec -> hsmmWorkingToConstrained() -> hsmmToHmm() ->
#' makeEmissionMatrix() -> scaledForward()` and returns the negative
#' log-likelihood, suitable for passing directly to [nlm()] (see
#' [fitHsmm()], which wraps exactly this call). Prefer
#' [safeHsmmNegativeLogLik()] as the actual optimizer objective: this
#' function raises an error on invalid parameter values, which aborts
#' optimization rather than steering the optimizer away from them.
#'
#' @param parVec Numeric working parameter vector.
#' @param x Numeric vector of observations (`NA` allowed for missing data).
#' @param modelTemplate A constrained HSMM model providing structural
#'   information not carried in `parVec` (see [hsmmWorkingToConstrained()]).
#' @param stationary Logical; whether the expanded-HMM stationary
#'   initialization (`TRUE`) or the aggregate-start initialization
#'   (`FALSE`) is used.
#' @param returnDetails Logical; if `TRUE`, returns the full intermediate
#'   pipeline instead of just the negative log-likelihood.
#'
#' @return If `returnDetails = FALSE` (the default), a single numeric value,
#'   the negative log-likelihood. If `returnDetails = TRUE`, a list
#'   containing `negativeLogLik`, `logLik`, `model` (the rebuilt constrained
#'   model), `hmmModel` (the expanded HMM), `emissionMatrix`, and
#'   `forwardResult`.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' workingStart <- hsmmConstrainedToWorking(testModel, stationary = TRUE)
#' nll <- hsmmNegativeLogLik(
#'   parVec = workingStart,
#'   x = rnorm(50),
#'   modelTemplate = testModel,
#'   stationary = TRUE
#' )
#'
#' @export
hsmmNegativeLogLik <- function(
    parVec,
    x,
    modelTemplate,
    stationary = TRUE,
    returnDetails = FALSE
) {
  model <- hsmmWorkingToConstrained(
    parVec = parVec,
    modelTemplate = modelTemplate,
    stationary = stationary
  )
  
  initialType <- if (stationary) {
    "stationaryHmm"
  } else {
    "aggregateStart"
  }
  
  hmmModel <- hsmmToHmm(
    model = model,
    initialType = initialType,
    checkRows = TRUE
  )
  
  emissionMatrix <- makeEmissionMatrix(
    x = x,
    stateParameters = hmmModel$stateParameters,
    emissionType = hmmModel$emissionType
  )
  
  forwardResult <- scaledForward(
    delta = hmmModel$delta,
    gamma = hmmModel$gamma,
    emissionMatrix = emissionMatrix,
    returnAlpha = returnDetails
  )
  
  if (returnDetails) {
    return(list(
      negativeLogLik = forwardResult$negativeLogLik,
      logLik = forwardResult$logLik,
      model = model,
      hmmModel = hmmModel,
      emissionMatrix = emissionMatrix,
      forwardResult = forwardResult
    ))
  }
  
  forwardResult$negativeLogLik
}

#' Safe negative log-likelihood wrapper for optimization
#'
#' Wraps [hsmmNegativeLogLik()] and returns a large finite penalty instead of
#' raising an error or returning a non-finite value, so that a numerical
#' optimizer exploring bad parameter regions is steered away from them
#' rather than aborted. This is the objective [fitHsmm()] actually passes to
#' [nlm()].
#'
#' @inheritParams hsmmNegativeLogLik
#' @param penalty Numeric penalty value returned on failure. When `NULL`
#'   (the default), it is `1e3 * length(x)`: clearly worse than any
#'   plausible negative log-likelihood, but not so large that `nlm()`'s
#'   finite-difference gradients explode when a penalized point sits next to
#'   a valid one. Scaling with the series length keeps it a few orders of
#'   magnitude above a typical fit rather than a hundred orders above it.
#'
#' @return A single numeric value: the negative log-likelihood, or
#'   `penalty` if evaluation failed or produced a non-finite value.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' workingStart <- hsmmConstrainedToWorking(testModel, stationary = TRUE)
#' nll <- safeHsmmNegativeLogLik(
#'   parVec = workingStart,
#'   x = rnorm(50),
#'   modelTemplate = testModel,
#'   stationary = TRUE
#' )
#'
#' @export
safeHsmmNegativeLogLik <- function(
    parVec,
    x,
    modelTemplate,
    stationary = TRUE,
    penalty = NULL
) {
  if (is.null(penalty)) {
    penalty <- 1e3 * max(1, length(x))
  }
  
  if (any(!is.finite(parVec))) {
    return(penalty)
  }
  
  value <- tryCatch(
    {
      nll <- hsmmNegativeLogLik(
        parVec = parVec,
        x = x,
        modelTemplate = modelTemplate,
        stationary = stationary
      )
      
      if (!is.finite(nll)) {
        penalty
      } else {
        nll
      }
    },
    error = function(errorCondition) {
      penalty
    }
  )
  
  value
}

# Internal helper for viterbi.R: log() that maps zero/invalid probabilities to
# -Inf instead of NaN, so Viterbi can represent impossible states/transitions
# exactly rather than through underflow.
safeLogProb <- function(probVec) {
  logProb <- rep(-Inf, length(probVec))
  validIndex <- probVec > 0 & is.finite(probVec)
  logProb[validIndex] <- log(probVec[validIndex])
  logProb
}
