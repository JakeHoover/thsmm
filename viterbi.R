#' Map expanded HMM substates back to HSMM aggregate states
#'
#' @param expandedStates Integer vector of expanded HMM state indices, e.g.
#'   `viterbiLog(...)$expandedStates`.
#' @param aggregateStarts,aggregateEnds Integer vectors giving the first and
#'   last expanded state in each HSMM aggregate (see [hsmmToHmm()]).
#'
#' @return An integer vector the same length as `expandedStates`, giving the
#'   corresponding HSMM state index for each.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' hmmModel <- hsmmToHmm(testModel)
#' expandedStatesToHsmmStates(
#'   expandedStates = c(1L, 2L, 11L),
#'   aggregateStarts = hmmModel$aggregateStarts,
#'   aggregateEnds = hmmModel$aggregateEnds
#' )
#'
#' @export
expandedStatesToHsmmStates <- function(expandedStates, aggregateStarts, aggregateEnds) {
  hsmmStates <- integer(length(expandedStates))
  
  for (stateIndex in seq_along(aggregateStarts)) {
    inAggregate <- expandedStates >= aggregateStarts[stateIndex] &
      expandedStates <= aggregateEnds[stateIndex]
    
    hsmmStates[inAggregate] <- stateIndex
  }
  
  hsmmStates
}


#' Core HMM-level Viterbi recursion
#'
#' Runs Viterbi decoding for the expanded HMM in log space. Impossible
#' transitions and impossible initial states are represented by `-Inf`, not
#' zero, so a genuinely unreachable state is distinguished from a merely
#' low-probability one; the backtrace raises an error rather than silently
#' returning a path that uses an impossible transition.
#'
#' This function has a true `O(nObs * nStates^2)` cost with two explicit
#' nested loops (time and state), unlike [scaledForward()], whose inner step
#' is a single matrix multiplication; it is consequently the most expensive
#' single computation in the package for long series with large state
#' aggregates.
#'
#' @param delta Numeric vector, the expanded HMM initial distribution.
#' @param gamma Numeric matrix, the expanded HMM transition matrix.
#' @param emissionMatrix Numeric matrix of state-dependent likelihoods, as
#'   returned by [makeEmissionMatrix()].
#'
#' @return A list containing `expandedStates` (the most likely expanded-HMM
#'   state path), `logProb` (log probability of the decoded path), `score`
#'   (the Viterbi score matrix), and `backPointer` (the backpointer matrix).
#'
#' @export
viterbiLog <- function(delta, gamma, emissionMatrix) {
  nObs <- nrow(emissionMatrix)
  nStates <- ncol(emissionMatrix)
  
  if (length(delta) != nStates) {
    stop("delta length does not match number of states.")
  }
  
  if (nrow(gamma) != nStates || ncol(gamma) != nStates) {
    stop("gamma dimensions do not match number of states.")
  }
  
  if (any(emissionMatrix < 0, na.rm = TRUE)) {
    stop("emissionMatrix contains negative values.")
  }
  
  logDelta <- safeLogProb(delta)
  
  logGamma <- matrix(-Inf, nrow = nStates, ncol = nStates)
  validGamma <- gamma > 0 & is.finite(gamma)
  logGamma[validGamma] <- log(gamma[validGamma])
  
  logEmission <- matrix(-Inf, nrow = nObs, ncol = nStates)
  validEmission <- emissionMatrix > 0 & is.finite(emissionMatrix)
  logEmission[validEmission] <- log(emissionMatrix[validEmission])
  
  score <- matrix(-Inf, nrow = nObs, ncol = nStates)
  backPointer <- matrix(NA_integer_, nrow = nObs, ncol = nStates)
  
  score[1, ] <- logDelta + logEmission[1, ]
  
  if (nObs >= 2L) {
    for (timeIndex in 2:nObs) {
      for (toState in seq_len(nStates)) {
        transitionScores <- score[timeIndex - 1, ] + logGamma[, toState]
        
        bestPreviousState <- which.max(transitionScores)
        bestScore <- transitionScores[bestPreviousState]
        
        score[timeIndex, toState] <- bestScore + logEmission[timeIndex, toState]
        
        # which.max() silently returns 1 when every candidate is -Inf, so record the
        # backpointer as missing instead of pointing at an unreachable predecessor.
        backPointer[timeIndex, toState] <- if (is.finite(bestScore)) {
          bestPreviousState
        } else {
          NA_integer_
        }
      }
    }
  }
  
  if (all(!is.finite(score[nObs, ]))) {
    stop(
      "Viterbi found no path with positive probability. ",
      "Check that delta, gamma and the emission parameters are compatible with the data."
    )
  }
  
  expandedStates <- integer(nObs)
  expandedStates[nObs] <- which.max(score[nObs, ])
  
  if (nObs >= 2L) {
    for (timeIndex in seq(from = nObs - 1L, to = 1L, by = -1L)) {
      previousState <- backPointer[
        timeIndex + 1L,
        expandedStates[timeIndex + 1L]
      ]
      
      if (is.na(previousState)) {
        stop("Viterbi backtrace reached an unreachable state at time ", timeIndex + 1L, ".")
      }
      
      expandedStates[timeIndex] <- previousState
    }
  }
  
  list(
    expandedStates = expandedStates,
    logProb = max(score[nObs, ]),
    score = score,
    backPointer = backPointer
  )
}


#' Decode the most likely HSMM state path
#'
#' Full HSMM-level Viterbi wrapper. Takes a **constrained model**, not a
#' working parameter vector — to decode with a fitted model, pass
#' `fit$fittedModel` from [fitHsmm()]; to decode from a working vector,
#' convert it first with [hsmmWorkingToConstrained()].
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()] or
#'   `fit$fittedModel` from [fitHsmm()].
#' @param x Numeric vector of observations.
#' @param initialType One of `"stationaryHmm"`, `"aggregateStart"`, or
#'   `"user"`. When `NULL` (the default), it is taken from `model$stationary`:
#'   `"stationaryHmm"` if the model is stationary, `"aggregateStart"`
#'   otherwise.
#' @param deltaHmm Optional numeric vector, a user-supplied expanded-HMM
#'   initial distribution. Required when `initialType = "user"`.
#' @param returnDetails Logical; if `TRUE`, also returns the constrained
#'   model, the expanded HMM, the emission matrix, and the raw
#'   [viterbiLog()] result.
#' @param checkRows Logical; passed through to [hsmmToHmm()].
#'
#' @return A list containing `states` (decoded HSMM state path),
#'   `expandedStates` (decoded expanded-HMM state path), and `logProb` (log
#'   probability of the decoded expanded-HMM path). With
#'   `returnDetails = TRUE`, also `model`, `hmmModel`, `emissionMatrix`, and
#'   `viterbiResult`.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 100, seed = 1)
#' vitResult <- hsmmViterbi(model = testModel, x = simData$x)
#'
#' @export
hsmmViterbi <- function(
    model,
    x,
    initialType = NULL,
    deltaHmm = NULL,
    returnDetails = FALSE,
    checkRows = TRUE
) {
  if (is.null(initialType)) {
    initialType <- if (isTRUE(model$stationary)) {
      "stationaryHmm"
    } else {
      "aggregateStart"
    }
  }
  
  validInitialTypes <- c(
    "stationaryHmm",
    "aggregateStart",
    "user"
  )
  
  if (!initialType %in% validInitialTypes) {
    stop(
      "initialType must be one of: ",
      paste(validInitialTypes, collapse = ", "),
      "."
    )
  }
  
  if (initialType == "aggregateStart" && is.null(model$delta)) {
    stop(
      "model$delta is required when initialType = 'aggregateStart'."
    )
  }
  
  if (initialType == "user" && is.null(deltaHmm)) {
    stop(
      "deltaHmm is required when initialType = 'user'."
    )
  }
  
  hmmModel <- hsmmToHmm(
    model = model,
    initialType = initialType,
    deltaHmm = deltaHmm,
    checkRows = checkRows
  )
  
  emissionMatrix <- makeEmissionMatrix(
    x = x,
    stateParameters = hmmModel$stateParameters,
    emissionType = hmmModel$emissionType,
    minDensity = 0
  )
  
  viterbiResult <- viterbiLog(
    delta = hmmModel$delta,
    gamma = hmmModel$gamma,
    emissionMatrix = emissionMatrix
  )
  
  hsmmStates <- expandedStatesToHsmmStates(
    expandedStates = viterbiResult$expandedStates,
    aggregateStarts = hmmModel$aggregateStarts,
    aggregateEnds = hmmModel$aggregateEnds
  )
  
  output <- list(
    states = hsmmStates,
    expandedStates = viterbiResult$expandedStates,
    logProb = viterbiResult$logProb
  )
  
  if (returnDetails) {
    output$model <- model
    output$hmmModel <- hmmModel
    output$emissionMatrix <- emissionMatrix
    output$viterbiResult <- viterbiResult
  }
  
  output
}


# ---- HSMM Fit Wrapper ----

