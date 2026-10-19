# Internal helpers: single-draw sampling and the state/emission simulation
# steps used by simulateHsmm() below.

sampleOne <- function(probVec) {
  sample.int(
    n = length(probVec),
    size = 1L,
    prob = probVec
  )
}


simulateMarkovChain <- function(delta, gamma, nObs) {
  nStates <- length(delta)
  states <- integer(nObs)
  
  states[1] <- sampleOne(delta)
  
  if (nObs >= 2L) {
    for (timeIndex in 2:nObs) {
      states[timeIndex] <- sampleOne(gamma[states[timeIndex - 1L], ])
    }
  }
  
  states
}


simulateEmissions <- function(expandedStates, stateParameters, emissionType) {
  nObs <- length(expandedStates)
  x <- numeric(nObs)
  
  if (emissionType == "normal") {
    mu <- stateParameters[, 1]
    sd <- stateParameters[, 2]
    
    if (any(sd <= 0) || any(!is.finite(sd))) {
      stop("Normal emission standard deviations must be positive and finite.")
    }
    
    for (timeIndex in seq_len(nObs)) {
      stateIndex <- expandedStates[timeIndex]
      x[timeIndex] <- rnorm(1L, mean = mu[stateIndex], sd = sd[stateIndex])
    }
  } else if (emissionType == "pois") {
    lambda <- stateParameters[, 1]
    
    if (any(lambda <= 0) || any(!is.finite(lambda))) {
      stop("Poisson emission rates must be positive and finite.")
    }
    
    for (timeIndex in seq_len(nObs)) {
      stateIndex <- expandedStates[timeIndex]
      x[timeIndex] <- rpois(1L, lambda = lambda[stateIndex])
    }
  } else {
    stop("Unsupported emissionType: ", emissionType)
  }
  
  x
}


#' Simulate data from an HSMM
#'
#' Generates synthetic data directly from a constructed HSMM model, using
#' the same expanded-HMM representation used for likelihood evaluation and
#' Viterbi decoding: HSMM model -> expanded HMM -> expanded state path ->
#' HSMM aggregate state path -> observations. This gives a convenient
#' end-to-end check of the full workflow (simulate -> evaluate likelihood
#' -> Viterbi decode) and is the standard way to generate test data for
#' recovery studies.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()].
#' @param nObs Integer, the number of observations to simulate.
#' @param initialType One of `"stationaryHmm"` (default), `"aggregateStart"`,
#'   or `"user"`; see [hsmmToHmm()].
#' @param deltaHmm Optional numeric vector, a user-supplied expanded-HMM
#'   initial distribution. Required when `initialType = "user"`.
#' @param seed Optional integer seed, set via [set.seed()] before
#'   simulating, for reproducibility.
#' @param returnDetails Logical; if `TRUE` (the default), also returns the
#'   original HSMM model and the expanded HMM object used for simulation.
#'
#' @return A list containing `x` (simulated observations: a vector for scalar
#'   emissions, an nObs x L matrix or nObs x L x d array of curves for
#'   functional ones), `states` (the
#'   HSMM-level state path), and `expandedStates` (the expanded-HMM state
#'   path). With `returnDetails = TRUE`, also `model` and `hmmModel`.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 500, seed = 123)
#' head(simData$x)
#' table(simData$states)
#'
#' @export
simulateHsmm <- function(
    model,
    nObs,
    initialType = c("stationaryHmm", "aggregateStart", "user"),
    deltaHmm = NULL,
    seed = NULL,
    returnDetails = TRUE
) {
  initialType <- match.arg(initialType)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  hmmModel <- hsmmToHmm(
    model = model,
    initialType = initialType,
    deltaHmm = deltaHmm,
    checkRows = TRUE
  )
  
  expandedStates <- simulateMarkovChain(
    delta = hmmModel$delta,
    gamma = hmmModel$gamma,
    nObs = nObs
  )
  
  states <- expandedStatesToHsmmStates(
    expandedStates = expandedStates,
    aggregateStarts = hmmModel$aggregateStarts,
    aggregateEnds = hmmModel$aggregateEnds
  )
  
  family <- getEmissionFamily(model$emissionType)

  if (family$isFunctional) {
    # Curves are generated per HSMM state from the unexpanded parameters;
    # emissionControl supplies the grid (or nPoints) and family settings.
    x <- family$simulate(states, model$emissionParameters, model$emissionControl)
  } else {
    x <- simulateEmissions(
      expandedStates = expandedStates,
      stateParameters = hmmModel$stateParameters,
      emissionType = hmmModel$emissionType
    )
  }
  
  output <- list(
    x = x,
    states = states,
    expandedStates = expandedStates
  )
  
  if (returnDetails) {
    output$model <- model
    output$hmmModel <- hmmModel
  }
  
  output
}


#' Realized dwell times from a simulated or decoded state path
#'
#' A simple diagnostic that reads realized dwell times off of a state path,
#' e.g. `simulateHsmm(...)$states` or `hsmmViterbi(...)$states`, by
#' run-length encoding it.
#'
#' @param states Integer (or otherwise atomic) vector of HSMM states.
#'
#' @return A data frame with columns `state` and `dwellTime`, one row per
#'   run.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 500, seed = 123)
#' stateRuns <- getStateRuns(simData$states)
#' \donttest{
#' boxplot(dwellTime ~ state, data = stateRuns)
#' }
#'
#' @export
getStateRuns <- function(states) {
  runInfo <- rle(states)
  
  data.frame(
    state = runInfo$values,
    dwellTime = runInfo$lengths
  )
}

