#' Build a constrained HSMM model object
#'
#' Constructs and validates a hidden semi-Markov model (HSMM) specification.
#' This is the entry point of the package: every other function in
#' \pkg{hsmmExpand} operates on the list returned here, either directly or
#' after it has been converted into an expanded HMM with [hsmmToHmm()].
#'
#' `buildHsmm()` performs basic dimension checks, validates and normalizes
#' `delta` whenever it is supplied, and infers `aggregateSizes` from
#' `sojournParameters` and `sojournType`. Use this constructor instead of
#' assembling the model list by hand.
#'
#' @param tpm A square numeric matrix, the embedded HSMM transition matrix.
#'   Rows must sum to one and the diagonal should be zero.
#' @param sojournParameters A list with one entry per HSMM state, giving the
#'   dwell-time parameters for that state. See [getDwellPmf()] for the
#'   parameter format expected for each `sojournType`.
#' @param sojournType A character vector with one entry per HSMM state,
#'   naming the dwell-time distribution family (e.g. `"geo"`, `"shiftp"`,
#'   `"nb"`, `"shiftnb"`, `"weibull"`, `"np"`, `"npg"`).
#' @param emissionParameters A numeric matrix with one row per HSMM state
#'   holding the state-dependent emission parameters.
#' @param emissionType A single string naming the emission distribution
#'   (`"normal"`, `"pois"`, or a functional family `"bmwd"`, `"ou"`,
#'   `"nonpar"`; functional families are fitted with [fitHsmmEM()] only).
#' @param delta Optional numeric vector, the HSMM-level initial distribution.
#'   Normalized whenever supplied. Required when `stationary = FALSE`;
#'   optional and ignored by `"stationaryHmm"` initialization otherwise.
#' @param stationary Logical; whether the model is intended to use
#'   stationary expanded-HMM initialization.
#' @param emissionControl Optional list of settings for the emission family
#'   (e.g. `grid`, `hurst` for `"bmwd"`; `norm`, `scale` for `"nonpar"`); see
#'   emissionFamilies.R. Not used by the scalar families.
#' @param fitShift Logical; whether shifted dwell-time offsets are to be
#'   estimated. Note that shift parameters cannot be optimized directly by
#'   [fitHsmm()] and must instead be profiled with [fitHsmmShiftGrid()].
#'
#' @return A list with class-free structure (a plain `list`) containing
#'   `nStates`, `tpm`, `delta`, `sojournParameters`, `sojournType`,
#'   `emissionParameters`, `emissionType`, `emissionControl`,
#'   `aggregateSizes`, `fitShift`, and
#'   `stationary`. `aggregateSizes` is inferred automatically and should not
#'   be supplied directly.
#'
#' @examples
#' # See exampleHsmmModel() for a runnable call to buildHsmm().
#' testModel <- exampleHsmmModel()
#' str(testModel)
#'
#' @export
buildHsmm <- function(
    tpm,
    sojournParameters,
    sojournType,
    emissionParameters,
    emissionType,
    delta = NULL,
    stationary = TRUE,
    fitShift = FALSE,
    emissionControl = NULL
) {
  nStates <- nrow(tpm)
  
  if (ncol(tpm) != nStates) {
    stop("tpm must be square.")
  }
  
  if (length(sojournParameters) != nStates) {
    stop("sojournParameters must have one entry per state.")
  }
  
  if (length(sojournType) != nStates) {
    stop("sojournType must have one entry per state.")
  }
  
  if (nrow(emissionParameters) != nStates) {
    stop("emissionParameters must have one row per state.")
  }
  
  aggregateSizes <- as.integer(mapply(
    FUN = getDwellSize,
    sojournParameters = sojournParameters,
    sojournType = sojournType
  ))
  
  if (!stationary && is.null(delta)) {
    stop("delta is required when stationary = FALSE.")
  }
  
  # Validate and normalise delta whenever it is supplied, not only when
  # stationary = FALSE, so a stored delta is always a proper distribution.
  if (!is.null(delta)) {
    if (length(delta) != nStates) {
      stop("delta must have length equal to the number of HSMM states.")
    }
    
    if (any(!is.finite(delta)) || any(delta < 0)) {
      stop("delta must contain finite non-negative values.")
    }
    
    if (sum(delta) <= 0) {
      stop("delta must have positive total mass.")
    }
    
    delta <- delta / sum(delta)
  }
  
  list(
    nStates = nStates,
    tpm = tpm,
    delta = delta,
    sojournParameters = sojournParameters,
    sojournType = sojournType,
    emissionParameters = emissionParameters,
    emissionType = emissionType,
    emissionControl = emissionControl,
    aggregateSizes = aggregateSizes,
    fitShift = fitShift,
    stationary = stationary
  )
}


# Internal helpers: aggregate index bookkeeping and expanded-HMM assembly.
# These are not exported; buildHsmm() and hsmmToHmm() are the public entry
# points that use them.

getAggregateStarts <- function(aggregateSizes) {
  # Base-R equivalent of c(1, head(cumsum(x) + 1, -1)). utils::head() would
  # pull in a dependency on utils for a single call, and leaving it
  # unimported trips codetools ("no visible global function definition").
  nAggregates <- length(aggregateSizes)
  c(1, (cumsum(aggregateSizes) + 1)[seq_len(nAggregates - 1L)])
}


getAggregateEnds <- function(aggregateSizes) {
  cumsum(aggregateSizes)
}


# tpm checks 
checkTransitionMatrix <- function(gamma, tol = 1e-8, name = "gamma") {
  if (!is.matrix(gamma)) {
    stop(name, " must be a matrix.")
  }
  
  if (nrow(gamma) != ncol(gamma)) {
    stop(name, " must be square.")
  }
  
  if (any(!is.finite(gamma))) {
    stop(name, " contains non-finite values.")
  }
  
  if (any(gamma < -tol)) {
    stop(name, " contains negative values below tolerance.")
  }
  
  rowSumError <- max(abs(rowSums(gamma) - 1))
  
  if (rowSumError > tol) {
    stop(
      name,
      " row sums are not equal to one. Maximum row-sum error: ",
      signif(rowSumError, 4)
    )
  }
  
  invisible(TRUE)
}


# build stationary distribution of HMM tpm 
getStationaryDistribution <- function(gamma, tol = 1e-10) {
  if (!is.matrix(gamma)) {
    stop("gamma must be a matrix.")
  }
  
  if (nrow(gamma) != ncol(gamma)) {
    stop("gamma must be square.")
  }
  
  if (any(!is.finite(gamma))) {
    stop("gamma contains non-finite values.")
  }
  
  gamma[gamma < 0] <- 0
  
  rowTotals <- rowSums(gamma)
  
  if (any(rowTotals <= 0)) {
    stop("gamma has at least one row with non-positive row sum.")
  }
  
  gamma <- gamma / rowTotals
  
  nStates <- nrow(gamma)
  
  # Solve delta (I - gamma + U) = 1, where U is a matrix of ones. This is the standard
  # rank-one correction and is better conditioned than replacing one balance equation
  # with the sum-to-one constraint.
  systemMatrix <- t(diag(nStates) - gamma + 1)
  systemVector <- rep(1, nStates)
  
  stationaryDistribution <- tryCatch(
    solve(systemMatrix, systemVector),
    error = function(errorCondition) {
      eigenFit <- eigen(t(gamma))
      eigenIndex <- which.min(abs(eigenFit$values - 1))
      candidate <- Re(eigenFit$vectors[, eigenIndex])
      
      candidate[candidate < 0] <- 0
      
      if (sum(candidate) <= tol) {
        candidate <- abs(Re(eigenFit$vectors[, eigenIndex]))
      }
      
      candidate / sum(candidate)
    }
  )
  
  # Clamp only genuinely negative values from round-off. Do not impose a positive
  # floor: tail substates of a large aggregate can have legitimate stationary
  # probabilities well below any fixed tolerance, and flooring them inflates them
  # by orders of magnitude and makes the objective non-smooth.
  stationaryDistribution <- pmax(stationaryDistribution, 0)
  
  totalMass <- sum(stationaryDistribution)
  
  if (!is.finite(totalMass) || totalMass <= 0) {
    stop("Stationary distribution could not be normalized.")
  }
  
  stationaryDistribution / totalMass
}


# expand emission distribution to state aggregates 
expandEmissionParameters <- function(emissionParameters, aggregateSizes, emissionType) {
  expandedParameters <- do.call(
    rbind,
    lapply(seq_along(aggregateSizes), function(stateIndex) {
      matrix(
        emissionParameters[stateIndex, ],
        nrow = aggregateSizes[stateIndex],
        ncol = ncol(emissionParameters),
        byrow = TRUE
      )
    })
  )
  
  if (emissionType == "normal") {
    colnames(expandedParameters) <- c("mean", "sd")
  }
  
  if (emissionType == "pois") {
    colnames(expandedParameters) <- "lambda"
  }
  
  expandedParameters
}

#' Convert an HSMM model into its expanded HMM representation
#'
#' Converts a constrained HSMM model (as built by [buildHsmm()]) into the
#' structured HMM obtained by replacing each HSMM state with a collection of
#' HMM substates ("aggregate"), with dwell-time hazards encoded in the
#' expanded transition matrix. This expanded representation is what the
#' likelihood ([scaledForward()], [hsmmNegativeLogLik()]), Viterbi decoding
#' ([hsmmViterbi()]), and simulation ([simulateHsmm()]) functions operate on.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()].
#' @param initialType One of `"stationaryHmm"` (default; the stationary
#'   distribution of the expanded transition matrix), `"aggregateStart"`
#'   (puts the HSMM-level initial distribution on the first substate of each
#'   aggregate), or `"user"` (uses the supplied `deltaHmm`).
#' @param deltaHmm Optional numeric vector, a user-supplied expanded-HMM
#'   initial distribution. Required when `initialType = "user"`.
#' @param checkRows Logical; whether to validate that `model$tpm` and the
#'   resulting expanded transition matrix are valid row-stochastic matrices.
#' @param tol Numeric tolerance used for the row-sum and diagonal checks.
#'
#' @return A list containing `gamma` (expanded HMM transition matrix),
#'   `stateParameters` (expanded state-dependent emission parameters),
#'   `emissionType`, `nHmmStates`, `delta` (expanded HMM initial
#'   distribution), `aggregateSizes`, `aggregateStarts`, `aggregateEnds`,
#'   `hazardValues`, `dwellPmfs`, and `initialType`.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' hmmModel <- hsmmToHmm(testModel, initialType = "stationaryHmm")
#' max(abs(rowSums(hmmModel$gamma) - 1))
#' sum(hmmModel$delta)
#'
#' @export
hsmmToHmm <- function(
    model,
    initialType = c("stationaryHmm", "aggregateStart", "user"),
    deltaHmm = NULL,
    checkRows = TRUE,
    tol = 1e-8
) {
  initialType <- match.arg(initialType)
  
  sojournType <- model$sojournType
  omega <- model$tpm
  dwellParameters <- model$sojournParameters
  nStates <- model$nStates
  emissionParameters <- model$emissionParameters
  emissionType <- model$emissionType
  deltaHsmm <- model$delta
  fitShift <- model$fitShift
  
  aggregateSizes <- as.integer(mapply(
    FUN = getDwellSize,
    sojournParameters = dwellParameters,
    sojournType = sojournType
  ))
  
  if (checkRows) {
    checkTransitionMatrix(omega, tol = tol, name = "omega")
    
    diagonalError <- max(abs(diag(omega)))
    
    if (diagonalError > tol) {
      stop(
        "omega should be the embedded HSMM transition matrix, so its diagonal should be zero. ",
        "Maximum diagonal entry: ",
        signif(diagonalError, 4)
      )
    }
  }
  
  totalStates <- sum(aggregateSizes)
  gamma <- matrix(0, nrow = totalStates, ncol = totalStates)
  
  aggregateStarts <- getAggregateStarts(aggregateSizes)
  aggregateEnds <- getAggregateEnds(aggregateSizes)
  
  dwellBlocks <- vector("list", nStates)
  hazardValues <- vector("list", nStates)
  dwellPmfs <- vector("list", nStates)
  
  for (stateIndex in seq_len(nStates)) {
    dwellBlock <- makeDwellBlock(
      sojournParameters = dwellParameters[[stateIndex]],
      sojournType = sojournType[stateIndex]
    )
    
    dwellBlocks[[stateIndex]] <- dwellBlock$block
    hazardValues[[stateIndex]] <- dwellBlock$hazard
    dwellPmfs[[stateIndex]] <- dwellBlock$dwellPmf
    
    rowIndex <- aggregateStarts[stateIndex]:aggregateEnds[stateIndex]
    colIndex <- aggregateStarts[stateIndex]:aggregateEnds[stateIndex]
    
    gamma[rowIndex, colIndex] <- dwellBlock$block
  }
  
  for (fromState in seq_len(nStates)) {
    rowIndex <- aggregateStarts[fromState]:aggregateEnds[fromState]
    
    for (toState in seq_len(nStates)) {
      if (fromState == toState) {
        next
      }
      
      toStart <- aggregateStarts[toState]
      gamma[rowIndex, toStart] <- omega[fromState, toState] * hazardValues[[fromState]]
    }
  }
  
  if (checkRows) {
    checkTransitionMatrix(gamma, tol = tol, name = "gamma")
  }
  
  stateParameters <- expandEmissionParameters(
    emissionParameters = emissionParameters,
    aggregateSizes = aggregateSizes,
    emissionType = emissionType
  )
  
  if (initialType == "aggregateStart") {
    if (is.null(deltaHsmm)) {
      stop("deltaHsmm is required when initialType = 'aggregateStart'.")
    }
    
    delta <- rep(0, totalStates)
    delta[aggregateStarts] <- deltaHsmm
  }
  
  if (initialType == "stationaryHmm") {
    delta <- getStationaryDistribution(gamma, tol = tol)
  }
  
  if (initialType == "user") {
    if (is.null(deltaHmm)) {
      stop("deltaHmm must be supplied when initialType = 'user'.")
    }
    
    if (length(deltaHmm) != totalStates) {
      stop("deltaHmm must have length equal to the number of expanded HMM states.")
    }
    
    if (any(deltaHmm < -tol)) {
      stop("deltaHmm contains negative values.")
    }
    
    delta <- pmax(deltaHmm, 0)
    delta <- delta / sum(delta)
  }
  
  list(
    gamma = gamma,
    stateParameters = stateParameters,
    emissionType = emissionType,
    nHmmStates = totalStates,
    delta = delta,
    fitShift = fitShift,
    aggregateSizes = aggregateSizes,
    aggregateStarts = aggregateStarts,
    aggregateEnds = aggregateEnds,
    hazardValues = hazardValues,
    dwellPmfs = dwellPmfs,
    initialType = initialType
  )
}
