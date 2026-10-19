# Internal parameter-transformation helpers.
#
# nlm() optimizes unconstrained vectors, so constrained model parameters
# (probabilities, simplices, positive scales) are mapped to and from an
# unconstrained "working" scale here. None of these are exported; the public
# entry points are makeIndexLayout(), hsmmConstrainedToWorking(), and
# hsmmWorkingToConstrained() below.

safeLog <- function(x, eps = 1e-12) {
  log(pmax(x, eps))
}

safeLogit <- function(x, eps = 1e-12) {
  x <- pmin(pmax(x, eps), 1 - eps)
  qlogis(x)
}

simplexToWorking <- function(probVec, eps = 1e-12) {
  probVec <- pmax(as.numeric(probVec), eps)
  probVec <- probVec / sum(probVec)
  log(probVec[-1] / probVec[1])
}

workingToSimplex <- function(workingVec) {
  rawVals <- exp(c(0, workingVec))
  rawVals / sum(rawVals)
}


tpmToWorking <- function(tpm, eps = 1e-12) {
  nStates <- nrow(tpm)
  
  if (nStates <= 2) {
    return(numeric(0))
  }
  
  workingList <- vector("list", nStates)
  
  for (stateIndex in seq_len(nStates)) {
    offDiagValues <- tpm[stateIndex, -stateIndex]
    workingList[[stateIndex]] <- simplexToWorking(offDiagValues, eps = eps)
  }
  
  unlist(workingList, use.names = FALSE)
}

workingToTpm <- function(workingVec, nStates) {
  if (nStates == 1) {
    return(matrix(1, 1, 1))
  }
  
  if (nStates == 2) {
    return(matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
  }
  
  expectedLength <- nStates * (nStates - 2)
  
  if (length(workingVec) != expectedLength) {
    stop(
      "workingVec has length ", length(workingVec),
      ", but expected ", expectedLength, " for the embedded TPM."
    )
  }
  
  tpm <- matrix(0, nrow = nStates, ncol = nStates)
  cursor <- 1L
  
  for (stateIndex in seq_len(nStates)) {
    rowWorking <- workingVec[cursor:(cursor + nStates - 3L)]
    cursor <- cursor + nStates - 2L
    
    offDiagValues <- workingToSimplex(rowWorking)
    tpm[stateIndex, -stateIndex] <- offDiagValues
  }
  
  tpm
}


getEmissionParamCount <- function(emissionType) {
  if (emissionType == "normal") {
    return(2L)
  }
  
  if (emissionType == "pois") {
    return(1L)
  }
  
  stop("Unsupported emissionType: ", emissionType)
}

emissionToWorking <- function(emissionParameters, emissionType) {
  emissionWorking <- emissionParameters
  
  if (emissionType == "normal") {
    emissionWorking[, 2] <- safeLog(emissionWorking[, 2])
  } else if (emissionType == "pois") {
    emissionWorking[, 1] <- safeLog(emissionWorking[, 1])
  } else {
    stop("Unsupported emissionType: ", emissionType)
  }
  
  as.numeric(t(emissionWorking))
}

workingToEmission <- function(workingVec, nStates, emissionType) {
  nEmissionParams <- getEmissionParamCount(emissionType)
  expectedLength <- nStates * nEmissionParams
  
  if (length(workingVec) != expectedLength) {
    stop(
      "workingVec has length ", length(workingVec),
      ", but expected ", expectedLength, " for emission parameters."
    )
  }
  
  emissionParameters <- matrix(
    workingVec,
    nrow = nStates,
    ncol = nEmissionParams,
    byrow = TRUE
  )
  
  if (emissionType == "normal") {
    emissionParameters[, 2] <- exp(emissionParameters[, 2])
    colnames(emissionParameters) <- c("mean", "sd")
  }
  
  if (emissionType == "pois") {
    emissionParameters[, 1] <- exp(emissionParameters[, 1])
    colnames(emissionParameters) <- "lambda"
  }
  
  emissionParameters
}


getSojournParamCount <- function(sojournParameters, sojournType, fitShift = FALSE) {
  if (sojournType == "npg") {
    return(length(sojournParameters))
  }
  
  if (sojournType == "np") {
    return(length(sojournParameters) - 1L)
  }
  
  if (sojournType == "geo") {
    return(1L)
  }
  
  if (sojournType == "pois") {
    return(1L)
  }
  
  if (sojournType == "nb") {
    return(2L)
  }
  
  if (sojournType == "shiftp") {
    return(ifelse(fitShift, 2L, 1L))
  }
  
  if (sojournType == "shiftnb") {
    return(ifelse(fitShift, 3L, 2L))
  }
  
  if (sojournType == "weibull") {
    return(2L)
  }
  
  stop("Unsupported sojournType: ", sojournType)
}

sojournToWorkingOne <- function(sojournParameters, sojournType, fitShift = FALSE) {
  if (sojournType == "npg") {
    headValues <- sojournParameters[-length(sojournParameters)]
    tailValue <- sojournParameters[length(sojournParameters)]
    
    # simplexToWorking() already drops the reference category, so headRaw of length
    # dwellSize yields dwellSize - 1 working values. Do not drop a second element.
    headRaw <- c(1 - sum(headValues), headValues)
    headWorking <- simplexToWorking(headRaw)
    tailWorking <- safeLogit(tailValue)
    
    return(c(headWorking, tailWorking))
  }
  
  if (sojournType == "np") {
    return(simplexToWorking(sojournParameters))
  }
  
  if (sojournType == "geo") {
    return(safeLogit(sojournParameters[1]))
  }
  
  if (sojournType == "pois") {
    return(safeLog(sojournParameters[1]))
  }
  
  if (sojournType == "nb") {
    sizeWorking <- safeLog(sojournParameters[1])
    probWorking <- safeLogit(sojournParameters[2])
    return(c(sizeWorking, probWorking))
  }
  
  if (sojournType == "shiftp") {
    lambdaWorking <- safeLog(sojournParameters[1])
    
    if (!fitShift) {
      return(lambdaWorking)
    }
    
    shiftValue <- sojournParameters[2]
    dwellSize <- sojournParameters[3]
    
    scaledShift <- (2 * shiftValue - (1 + dwellSize)) / (dwellSize - 1)
    scaledShift <- pmin(pmax(scaledShift, -1 + 1e-12), 1 - 1e-12)
    
    return(c(lambdaWorking, atanh(scaledShift)))
  }
  
  if (sojournType == "shiftnb") {
    sizeWorking <- safeLog(sojournParameters[1])
    probWorking <- safeLogit(sojournParameters[2])
    
    if (!fitShift) {
      return(c(sizeWorking, probWorking))
    }
    
    shiftValue <- sojournParameters[3]
    dwellSize <- sojournParameters[4]
    
    scaledShift <- (2 * shiftValue - (1 + dwellSize)) / (dwellSize - 1)
    scaledShift <- pmin(pmax(scaledShift, -1 + 1e-12), 1 - 1e-12)
    
    return(c(sizeWorking, probWorking, atanh(scaledShift)))
  }
  
  if (sojournType == "weibull") {
    shapeWorking <- safeLog(sojournParameters[1])
    scaleWorking <- safeLog(sojournParameters[2])
    return(c(shapeWorking, scaleWorking))
  }
  
  stop("Unsupported sojournType: ", sojournType)
}

workingToSojournOne <- function(workingVec, sojournTemplate, sojournType, fitShift = FALSE) {
  if (sojournType == "npg") {
    headWorking <- workingVec[-length(workingVec)]
    tailWorking <- workingVec[length(workingVec)]
    
    headRaw <- exp(headWorking)
    headValues <- headRaw / (1 + sum(headRaw))
    tailValue <- plogis(tailWorking)
    
    return(c(headValues, tailValue))
  }
  
  if (sojournType == "np") {
    return(workingToSimplex(workingVec))
  }
  
  if (sojournType == "geo") {
    return(plogis(workingVec[1]))
  }
  
  if (sojournType == "pois") {
    lambdaValue <- exp(workingVec[1])
    dwellSize <- sojournTemplate[2]
    return(c(lambdaValue, dwellSize))
  }
  
  if (sojournType == "nb") {
    sizeValue <- exp(workingVec[1])
    probValue <- plogis(workingVec[2])
    dwellSize <- sojournTemplate[3]
    return(c(sizeValue, probValue, dwellSize))
  }
  
  if (sojournType == "shiftp") {
    lambdaValue <- exp(workingVec[1])
    dwellSize <- sojournTemplate[3]
    
    if (fitShift) {
      shiftValue <- tanh(workingVec[2])
      shiftValue <- 0.5 * ((dwellSize - 1) * shiftValue + (1 + dwellSize))
    } else {
      shiftValue <- sojournTemplate[2]
    }
    
    return(c(lambdaValue, shiftValue, dwellSize))
  }
  
  if (sojournType == "shiftnb") {
    sizeValue <- exp(workingVec[1])
    probValue <- plogis(workingVec[2])
    dwellSize <- sojournTemplate[4]
    
    if (fitShift) {
      shiftValue <- tanh(workingVec[3])
      shiftValue <- 0.5 * ((dwellSize - 1) * shiftValue + (1 + dwellSize))
    } else {
      shiftValue <- sojournTemplate[3]
    }
    
    return(c(sizeValue, probValue, shiftValue, dwellSize))
  }
  
  if (sojournType == "weibull") {
    shapeValue <- exp(workingVec[1])
    scaleValue <- exp(workingVec[2])
    dwellSize <- sojournTemplate[3]
    
    return(c(shape = shapeValue, scale = scaleValue, dwellSize = dwellSize))
  }
  
  stop("Unsupported sojournType: ", sojournType)
}


#' Build index bookkeeping for the working parameter vector
#'
#' Computes explicit index ranges for each section of the unconstrained
#' working parameter vector used by [hsmmConstrainedToWorking()],
#' [hsmmWorkingToConstrained()], and the optimizer-facing likelihood
#' functions. The working vector is ordered as transition parameters,
#' emission parameters, sojourn parameters, and (when `stationary = FALSE`)
#' initial-distribution parameters.
#'
#' Also useful on its own for computing the number of free parameters in a
#' model, e.g. for AIC/BIC model comparison — different sojourn families
#' have very different parameter counts, so raw negative log-likelihoods are
#' not comparable across models without this.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()].
#' @param stationary Logical; whether the initial-distribution block is
#'   included in the working vector.
#'
#' @return A list containing `counts` (named integer vector of section
#'   sizes), `sections` (named list of index vectors into the working
#'   vector), `sojournCounts` (parameter count per HSMM state), and
#'   `totalCount` (the total working-vector length).
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' layout <- makeIndexLayout(testModel, stationary = TRUE)
#' layout$totalCount
#'
#' @export
makeIndexLayout <- function(model, stationary = TRUE) {
  nStates <- model$nStates
  sojournTypes <- model$sojournType
  dwellParameters <- model$sojournParameters
  fitShift <- model$fitShift
  emissionType <- model$emissionType
  
  tpmCount <- ifelse(nStates > 2, nStates * (nStates - 2), 0L)
  emissionCount <- nStates * getEmissionParamCount(emissionType)
  
  sojournCounts <- mapply(
    FUN = getSojournParamCount,
    sojournParameters = dwellParameters,
    sojournType = sojournTypes,
    MoreArgs = list(fitShift = fitShift)
  )
  
  sojournCount <- sum(sojournCounts)
  deltaCount <- ifelse(stationary, 0L, nStates - 1L)
  
  counts <- c(
    tpm = tpmCount,
    emission = emissionCount,
    sojourn = sojournCount,
    delta = deltaCount
  )
  
  sections <- list()
  cursor <- 1L
  
  for (sectionName in names(counts)) {
    sectionCount <- as.integer(counts[sectionName])
    
    if (sectionCount == 0L) {
      sections[[sectionName]] <- integer(0)
    } else {
      sections[[sectionName]] <- cursor:(cursor + sectionCount - 1L)
      cursor <- cursor + sectionCount
    }
  }
  
  list(
    counts = as.integer(counts),
    sections = sections,
    sojournCounts = as.integer(sojournCounts),
    totalCount = sum(counts)
  )
}


# Splits a flat vector into a list of sub-vectors given a vector of section
# lengths. Used to pull the per-state sojourn working parameters back apart.
splitByCounts <- function(x, counts) {
  output <- vector("list", length(counts))
  cursor <- 1L
  
  for (index in seq_along(counts)) {
    count <- counts[index]
    
    if (count == 0L) {
      output[[index]] <- numeric(0)
    } else {
      output[[index]] <- x[cursor:(cursor + count - 1L)]
      cursor <- cursor + count
    }
  }
  
  output
}

#' Convert a constrained HSMM model to an unconstrained working vector
#'
#' Flattens a constrained HSMM model into the unconstrained working
#' parameter vector expected by [hsmmNegativeLogLik()] / [nlm()], via
#' [makeIndexLayout()].
#'
#' When `stationary = TRUE`, the initial distribution is not included in the
#' returned vector; the expanded-HMM stationary distribution is instead
#' computed from the expanded transition matrix each time the model is
#' rebuilt from a working vector. When `stationary = FALSE`, the HSMM-level
#' `delta` is included as simplex working parameters.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()].
#' @param stationary Logical; must match the `stationary` value used
#'   everywhere else in the fitting pipeline for this model.
#'
#' @return A numeric vector, the unconstrained working parameters.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' workingStart <- hsmmConstrainedToWorking(testModel, stationary = TRUE)
#'
#' @export
hsmmConstrainedToWorking <- function(model, stationary = TRUE) {
  layout <- makeIndexLayout(model, stationary = stationary)
  
  tpmWorking <- tpmToWorking(model$tpm)
  
  emissionWorking <- emissionToWorking(
    emissionParameters = model$emissionParameters,
    emissionType = model$emissionType
  )
  
  sojournWorkingList <- mapply(
    FUN = sojournToWorkingOne,
    sojournParameters = model$sojournParameters,
    sojournType = model$sojournType,
    MoreArgs = list(fitShift = model$fitShift),
    SIMPLIFY = FALSE
  )
  
  sojournWorking <- unlist(sojournWorkingList, use.names = FALSE)
  
  if (stationary) {
    deltaWorking <- numeric(0)
  } else {
    if (is.null(model$delta)) {
      stop("model$delta is required when stationary = FALSE.")
    }
    
    deltaWorking <- simplexToWorking(model$delta)
  }
  
  parVec <- c(
    tpmWorking,
    emissionWorking,
    sojournWorking,
    deltaWorking
  )
  
  if (length(parVec) != layout$totalCount) {
    stop(
      "Parameter vector has length ", length(parVec),
      ", but expected ", layout$totalCount, "."
    )
  }
  
  parVec
}


#' Convert a working parameter vector back into a constrained HSMM model
#'
#' The inverse of [hsmmConstrainedToWorking()]: rebuilds a full constrained
#' HSMM model from an unconstrained working vector, using `modelTemplate`
#' for structural information that is not estimated (number of states,
#' dwell-time families, fixed dwell-size truncation levels, fixed shift
#' values, and emission distribution type).
#'
#' @param parVec Numeric working parameter vector, as produced by
#'   [hsmmConstrainedToWorking()] or an optimizer such as [nlm()].
#' @param modelTemplate A constrained HSMM model providing structural
#'   information not carried in `parVec`.
#' @param stationary Logical; must match the `stationary` value `parVec`
#'   was built or optimized under.
#'
#' @return A constrained HSMM model (see [buildHsmm()]).
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' workingStart <- hsmmConstrainedToWorking(testModel, stationary = TRUE)
#' rebuiltModel <- hsmmWorkingToConstrained(workingStart, testModel, stationary = TRUE)
#'
#' @export
hsmmWorkingToConstrained <- function(parVec, modelTemplate, stationary = TRUE) {
  layout <- makeIndexLayout(modelTemplate, stationary = stationary)
  
  if (length(parVec) != layout$totalCount) {
    stop(
      "parVec has length ", length(parVec),
      ", but expected ", layout$totalCount, "."
    )
  }
  
  sections <- layout$sections
  nStates <- modelTemplate$nStates
  
  tpm <- workingToTpm(
    workingVec = parVec[sections$tpm],
    nStates = nStates
  )
  
  emissionParameters <- workingToEmission(
    workingVec = parVec[sections$emission],
    nStates = nStates,
    emissionType = modelTemplate$emissionType
  )
  
  sojournWorkingList <- splitByCounts(
    x = parVec[sections$sojourn],
    counts = layout$sojournCounts
  )
  
  sojournParameters <- vector("list", nStates)
  
  for (stateIndex in seq_len(nStates)) {
    sojournParameters[[stateIndex]] <- workingToSojournOne(
      workingVec = sojournWorkingList[[stateIndex]],
      sojournTemplate = modelTemplate$sojournParameters[[stateIndex]],
      sojournType = modelTemplate$sojournType[stateIndex],
      fitShift = modelTemplate$fitShift
    )
  }
  
  if (stationary) {
    delta <- NULL
  } else {
    delta <- workingToSimplex(parVec[sections$delta])
  }
  
  buildHsmm(
    tpm = tpm,
    sojournParameters = sojournParameters,
    sojournType = modelTemplate$sojournType,
    emissionParameters = emissionParameters,
    emissionType = modelTemplate$emissionType,
    delta = delta,
    stationary = stationary,
    fitShift = modelTemplate$fitShift
  )
}
