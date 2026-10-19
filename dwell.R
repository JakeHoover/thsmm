#' Infer the number of expanded HMM substates for a dwell-time specification
#'
#' Determines the size of the state aggregate (the number of HMM substates)
#' needed to represent one HSMM state's dwell-time distribution, based on its
#' `sojournType` and `sojournParameters`. Used internally by [buildHsmm()]
#' and [hsmmToHmm()]; also useful on its own when constructing model
#' specifications by hand.
#'
#' @param sojournParameters Numeric vector of dwell-time parameters for one
#'   state; the expected layout depends on `sojournType` (see [getDwellPmf()]).
#' @param sojournType A single string naming the dwell-time family: one of
#'   `"geo"`, `"pois"`, `"shiftp"`, `"nb"`, `"shiftnb"`, `"weibull"`, `"np"`,
#'   or `"npg"`.
#'
#' @return A single integer, the number of expanded HMM substates.
#'
#' @examples
#' getDwellSize(c(lambda = 4, shift = 1, dwellSize = 10), "shiftp")
#' getDwellSize(c(0.1, 0.4, 0.3, 0.2), "np")
#' getDwellSize(c(prob = 0.25), "geo")
#'
#' @export
getDwellSize <- function(sojournParameters, sojournType) {
  if (sojournType == "pois") {
    return(as.integer(sojournParameters[2]))
  }
  
  if (sojournType == "shiftp") {
    return(as.integer(sojournParameters[3]))
  }
  
  if (sojournType == "nb") {
    return(as.integer(sojournParameters[3]))
  }
  
  if (sojournType == "shiftnb") {
    return(as.integer(sojournParameters[4]))
  }
  
  if (sojournType == "geo") {
    return(1L)
  }
  
  if (sojournType == "weibull") {
    return(as.integer(sojournParameters[3]))
  }
  
  if (sojournType %in% c("np", "npg")) {
    return(length(sojournParameters))
  }
  
  stop("Unsupported sojournType: ", sojournType)
}


#' Evaluate a dwell-time probability mass function
#'
#' Evaluates the dwell-time pmf over the truncated support used by the
#' expanded HMM, for one of the supported `sojournType` families.
#'
#' @param sojournParameters Numeric vector of dwell-time parameters, in the
#'   format expected for `sojournType`:
#'   \describe{
#'     \item{`"geo"`}{`c(prob)`}
#'     \item{`"pois"`}{`c(lambda, dwellSize)`}
#'     \item{`"shiftp"`}{`c(lambda, shift, dwellSize)`}
#'     \item{`"nb"`}{`c(size, prob, dwellSize)`}
#'     \item{`"shiftnb"`}{`c(size, prob, shift, dwellSize)`}
#'     \item{`"weibull"`}{`c(shape, scale, dwellSize)`}
#'     \item{`"np"`}{a probability vector summing to one}
#'     \item{`"npg"`}{`c(d(1), ..., d(m - 1), tailHazard)`, where the head
#'       entries are unstructured probabilities summing to less than one and
#'       the final entry is the geometric tail *hazard*, not a probability}
#'   }
#' @param sojournType A single string naming the dwell-time family.
#'
#' @return A numeric vector, the dwell-time pmf over its truncated support.
#'   For families with infinite support, this is the *nominal* pmf on
#'   `1:dwellSize`; the distribution the expanded HMM actually encodes also
#'   has a geometric tail beyond `dwellSize` and is given by
#'   [getExpandedDwellPmf()] instead.
#'
#' @export
getDwellPmf <- function(sojournParameters, sojournType) {
  dwellSize <- getDwellSize(sojournParameters, sojournType)
  dwellTimes <- seq_len(dwellSize)
  
  if (sojournType == "np") {
    dwellPmf <- as.numeric(sojournParameters)
  }
  
  else if (sojournType == "npg") {
    # The first dwellSize - 1 entries are unstructured probabilities d(1), ..., d(m - 1).
    # The last entry is the hazard of the geometric tail, not a probability mass.
    # The mass placed on the tail state is therefore (1 - sum(head)) * tail, which makes
    # dwellPmfToHazard() recover hazard[m] = tail exactly.
    headValues <- as.numeric(sojournParameters[-dwellSize])
    tailValue <- as.numeric(sojournParameters[dwellSize])
    
    headMass <- sum(headValues)
    
    if (headMass >= 1) {
      stop("npg head probabilities must sum to less than one.")
    }
    
    dwellPmf <- c(headValues, (1 - headMass) * tailValue)
  }
  
  else if (sojournType == "geo") {
    dwellPmf <- as.numeric(sojournParameters[1])
  }
  
  else if (sojournType == "pois") {
    lambda <- sojournParameters[1]
    dwellPmf <- dpois(dwellTimes - 1, lambda = lambda)
  }
  
  else if (sojournType == "nb") {
    size <- sojournParameters[1]
    prob <- sojournParameters[2]
    dwellPmf <- dnbinom(dwellTimes - 1, size = size, prob = prob)
  }
  
  else if (sojournType == "shiftp") {
    lambda <- sojournParameters[1]
    shift <- floor(sojournParameters[2])
    
    if (shift >= dwellSize) {
      warning(
        "shift = ", shift,
        " is greater than or equal to dwellSize = ", dwellSize,
        ". Setting shift to dwellSize - 1."
      )
      shift <- dwellSize - 1
    }
    
    dwellPmf <- ifelse(
      dwellTimes < shift,
      0,
      dpois(dwellTimes - shift, lambda = lambda)
    )
  }
  
  else if (sojournType == "shiftnb") {
    size <- sojournParameters[1]
    prob <- sojournParameters[2]
    shift <- floor(sojournParameters[3])
    
    if (shift >= dwellSize) {
      warning(
        "shift = ", shift,
        " is greater than or equal to dwellSize = ", dwellSize,
        ". Setting shift to dwellSize - 1."
      )
      shift <- dwellSize - 1
    }
    
    dwellPmf <- ifelse(
      dwellTimes < shift,
      0,
      dnbinom(dwellTimes - shift, size = size, prob = prob)
    )
  }
  
  else if (sojournType == "weibull") {
    shape <- sojournParameters[1]
    scale <- sojournParameters[2]
    
    dwellPmf <- pweibull(
      dwellTimes,
      shape = shape,
      scale = scale
    ) - pweibull(
      dwellTimes - 1,
      shape = shape,
      scale = scale
    )
  }
  
  else {
    stop("Unsupported sojournType: ", sojournType)
  }
  
  as.numeric(dwellPmf)
}


#' Dwell-time survival function
#'
#' Returns the dwell-time survival function `S(r) = 1 - F(r - 1) = P(D >= r)`
#' over the truncated support, computed from the upper tail of the
#' parametric CDF (e.g. `ppois(..., lower.tail = FALSE)`) rather than by
#' subtracting a cumulative sum. Subtraction loses all precision once the
#' survival probability is small, which produces hazards above one and then
#' `NaN` for aggregates that reach far into the tail.
#'
#' @inheritParams getDwellPmf
#'
#' @return A numeric vector of survival probabilities, or `NULL` for
#'   `"np"` and `"npg"`, where the support is short and the cumulative sum
#'   used by [dwellPmfToHazard()] is exact.
#'
#' @export
getDwellSurvival <- function(sojournParameters, sojournType) {
  dwellSize <- getDwellSize(sojournParameters, sojournType)
  dwellTimes <- seq_len(dwellSize)
  
  if (sojournType %in% c("np", "npg")) {
    return(NULL)
  }
  
  if (sojournType == "geo") {
    return(rep(1, dwellSize))
  }
  
  if (sojournType == "pois") {
    lambda <- sojournParameters[1]
    return(ppois(dwellTimes - 2, lambda = lambda, lower.tail = FALSE))
  }
  
  if (sojournType == "shiftp") {
    lambda <- sojournParameters[1]
    shift <- floor(sojournParameters[2])
    shift <- min(shift, dwellSize - 1)
    return(ppois(dwellTimes - shift - 1, lambda = lambda, lower.tail = FALSE))
  }
  
  if (sojournType == "nb") {
    size <- sojournParameters[1]
    prob <- sojournParameters[2]
    return(pnbinom(dwellTimes - 2, size = size, prob = prob, lower.tail = FALSE))
  }
  
  if (sojournType == "shiftnb") {
    size <- sojournParameters[1]
    prob <- sojournParameters[2]
    shift <- floor(sojournParameters[3])
    shift <- min(shift, dwellSize - 1)
    return(pnbinom(
      dwellTimes - shift - 1,
      size = size,
      prob = prob,
      lower.tail = FALSE
    ))
  }
  
  if (sojournType == "weibull") {
    shape <- sojournParameters[1]
    scale <- sojournParameters[2]
    return(pweibull(
      dwellTimes - 1,
      shape = shape,
      scale = scale,
      lower.tail = FALSE
    ))
  }
  
  stop("Unsupported sojournType: ", sojournType)
}


#' Convert a dwell-time pmf into hazard values
#'
#' Converts a dwell-time probability vector into discrete hazard values
#' `c_k(r) = p_k(r) / (1 - F_k(r - 1))`, which is how dwell-time
#' distributions are actually encoded in the expanded HMM transition matrix.
#'
#' @param dwellPmf Numeric vector, the dwell-time pmf over its support.
#' @param forceFiniteSupport Logical; if `TRUE`, the last hazard is set to
#'   one. Appropriate for finite-support dwell-time distributions such as
#'   `"np"`.
#' @param eps Numeric tolerance for treating small survival probabilities as
#'   zero.
#' @param survivalBefore Optional numeric vector of survival probabilities.
#'   If supplied (e.g. from [getDwellSurvival()]), it is used directly
#'   instead of being reconstructed from a cumulative sum of `dwellPmf`, which
#'   is the numerically stable choice for heavy-tailed families.
#'
#' @return A numeric vector of hazard values in `[0, 1]`, the same length as
#'   `dwellPmf`.
#'
#' @export
dwellPmfToHazard <- function(dwellPmf,
                             forceFiniteSupport = FALSE,
                             eps = 1e-12,
                             survivalBefore = NULL) {
  dwellPmf <- as.numeric(dwellPmf)
  dwellSize <- length(dwellPmf)
  
  if (any(!is.finite(dwellPmf))) {
    stop("dwellPmf contains non-finite values.")
  }
  
  if (any(dwellPmf < -eps)) {
    stop("dwellPmf contains negative values.")
  }
  
  dwellPmf <- pmax(dwellPmf, 0)
  
  if (is.null(survivalBefore)) {
    survivalBefore <- 1 - c(0, cumsum(dwellPmf[-dwellSize]))
  } else {
    survivalBefore <- as.numeric(survivalBefore)
    
    if (length(survivalBefore) != dwellSize) {
      stop("survivalBefore must have the same length as dwellPmf.")
    }
    
    if (any(!is.finite(survivalBefore))) {
      stop("survivalBefore contains non-finite values.")
    }
  }
  
  hazard <- ifelse(
    survivalBefore > eps,
    dwellPmf / survivalBefore,
    1
  )
  
  hazard <- pmin(pmax(hazard, 0), 1)
  
  if (forceFiniteSupport) {
    hazard[dwellSize] <- 1
  }
  
  hazard
}


#' Build the within-aggregate transition block for one HSMM state
#'
#' Creates the diagonal block of the expanded HMM transition matrix that
#' governs remaining inside one HSMM state's aggregate. The missing row mass
#' is assigned to transitions out of the aggregate when the full expanded
#' transition matrix is assembled by [hsmmToHmm()].
#'
#' @inheritParams dwellPmfToHazard
#' @param sojournParameters,sojournType Dwell-time specification; see
#'   [getDwellPmf()] for the expected parameter layout.
#'
#' @return A list containing `block` (the within-aggregate transition
#'   block), `hazard` (the dwell-time hazard values), and `dwellPmf` (the
#'   dwell-time probability values over the aggregate support).
#'
#' @export
makeDwellBlock <- function(sojournParameters, sojournType, eps = 1e-12) {
  dwellPmf <- getDwellPmf(sojournParameters, sojournType)
  dwellSize <- length(dwellPmf)
  
  survivalBefore <- getDwellSurvival(sojournParameters, sojournType)
  
  forceFiniteSupport <- sojournType == "np"
  hazard <- dwellPmfToHazard(
    dwellPmf = dwellPmf,
    forceFiniteSupport = forceFiniteSupport,
    eps = eps,
    survivalBefore = survivalBefore
  )
  
  block <- matrix(0, nrow = dwellSize, ncol = dwellSize)
  
  if (dwellSize == 1) {
    block[1, 1] <- 1 - hazard[1]
  } else {
    block[cbind(seq_len(dwellSize - 1), 2:dwellSize)] <- 1 - hazard[-dwellSize]
    block[dwellSize, dwellSize] <- 1 - hazard[dwellSize]
  }
  
  list(
    block = block,
    hazard = hazard,
    dwellPmf = dwellPmf
  )
}


#' Dwell-time distribution actually represented by the expanded HMM
#'
#' Returns `p*(r)`, the dwell-time distribution from Langrock and Zucchini
#' (2011), Proposition 2: it equals the nominal dwell-time pmf for
#' `r <= dwellSize` and decays geometrically beyond it. This is the pmf the
#' expanded HMM actually encodes, and the correct object to use for plotting
#' fitted dwell distributions and for computing distances between them,
#' because the nominal untruncated parametric pmf disagrees with it whenever
#' `dwellSize` truncates a distribution with meaningful tail mass.
#'
#' Because this works from the hazard values alone, it covers every sojourn
#' family with no per-family branching, and stays correct as new families
#' are added.
#'
#' @param hazard Numeric vector of hazard values, e.g. `makeDwellBlock(...)$hazard`.
#' @param maxDwell Integer, how far out to evaluate the pmf. Defaults to
#'   `length(hazard)` (no geometric tail).
#'
#' @return A numeric vector of length `maxDwell`. The sum should be close to
#'   one for a well-chosen `maxDwell`.
#'
#' @examples
#' dwellBlock <- makeDwellBlock(c(lambda = 4, shift = 1, dwellSize = 10), "shiftp")
#' expandedPmf <- getExpandedDwellPmf(hazard = dwellBlock$hazard, maxDwell = 100)
#' sum(expandedPmf)
#'
#' @export
getExpandedDwellPmf <- function(hazard, maxDwell = NULL) {
  hazard <- as.numeric(hazard)
  dwellSize <- length(hazard)
  
  if (dwellSize < 1L) {
    stop("hazard must have at least one element.")
  }
  
  if (any(!is.finite(hazard)) || any(hazard < 0) || any(hazard > 1)) {
    stop("hazard values must be finite and lie in [0, 1].")
  }
  
  if (is.null(maxDwell)) {
    maxDwell <- dwellSize
  }
  
  maxDwell <- as.integer(maxDwell)
  
  if (maxDwell < 1L) {
    stop("maxDwell must be a positive integer.")
  }
  
  # survivalBefore[r] = P(D >= r) = prod_{s < r} (1 - hazard[s])
  survivalBefore <- c(1, cumprod(1 - hazard))
  
  dwellPmf <- numeric(maxDwell)
  headIndex <- seq_len(min(dwellSize, maxDwell))
  dwellPmf[headIndex] <- hazard[headIndex] * survivalBefore[headIndex]
  
  if (maxDwell > dwellSize) {
    tailIndex <- (dwellSize + 1L):maxDwell
    dwellPmf[tailIndex] <- dwellPmf[dwellSize] *
      (1 - hazard[dwellSize])^(tailIndex - dwellSize)
  }
  
  dwellPmf
}


#' Expanded dwell-time pmf for one state of a fitted or constructed model
#'
#' Convenience wrapper around [makeDwellBlock()] and [getExpandedDwellPmf()]
#' for one state of an HSMM model. Use this for plotting fitted dwell
#' distributions and for computing distances between them; re-evaluating the
#' nominal parametric pmf instead (via [getDwellPmf()]) will disagree with
#' the fitted model wherever `dwellSize` truncates a heavy tail.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()] or
#'   `fit$fittedModel` from [fitHsmm()].
#' @param stateIndex Integer index of the HSMM state to evaluate.
#' @param maxDwell Integer, how far out to evaluate the pmf. Defaults to the
#'   state's own aggregate size (no geometric tail); pass a larger value to
#'   capture tail mass beyond the truncation.
#'
#' @return A numeric vector of length `maxDwell`, the dwell-time pmf.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' dwellPmf <- getModelDwellPmf(testModel, stateIndex = 1, maxDwell = 100)
#'
#' @export
getModelDwellPmf <- function(model, stateIndex, maxDwell = NULL) {
  dwellBlock <- makeDwellBlock(
    sojournParameters = model$sojournParameters[[stateIndex]],
    sojournType = model$sojournType[stateIndex]
  )
  
  getExpandedDwellPmf(
    hazard = dwellBlock$hazard,
    maxDwell = maxDwell
  )
}

