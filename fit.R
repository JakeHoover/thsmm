#' Fit an HSMM by maximum likelihood
#'
#' Runs [nlm()] on the unconstrained working parameter vector and returns
#' the fitted model in both working and constrained form. This is the
#' primary fitting entry point of the package.
#'
#' When `initialType` is `NULL` it is derived from `stationary`, so the
#' reported `fittedHmm` always uses the same initial distribution the
#' likelihood was evaluated under. A warning is raised when `nlm()` returns
#' convergence code 3, 4, or 5; code 4 means the iteration limit was
#' reached — raise `iterlim`, adjust `stepmax`, or try a different start.
#'
#' `fitShift = TRUE` on `model` does **not** work with this function: the
#' dwell-time pmf uses `floor(shift)`, so the likelihood is a step function
#' of the shift working parameter, its numerical gradient is zero almost
#' everywhere, and `nlm()` cannot move it. Use [fitHsmmShiftGrid()] to
#' profile over integer shifts instead; `fitHsmm()` warns when it sees
#' `fitShift = TRUE`.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()],
#'   giving the starting values and structural template for the fit.
#' @param dat Numeric vector of observations (`NA` allowed for missing
#'   data).
#' @param stationary Logical; whether to use the expanded-HMM stationary
#'   initialization (`TRUE`) or the aggregate-start initialization
#'   (`FALSE`) for likelihood evaluation.
#' @param initialType Initial-distribution mode used to build the reported
#'   `fittedHmm`; see [hsmmToHmm()]. Derived from `stationary` when `NULL`.
#' @param workingStart Optional numeric starting vector on the working
#'   scale. Computed from `model` via [hsmmConstrainedToWorking()] when
#'   `NULL`.
#' @param safeLikelihood Logical; whether to optimize
#'   [safeHsmmNegativeLogLik()] (recommended) instead of the raw
#'   [hsmmNegativeLogLik()].
#' @param iterlim,stepmax Passed to [nlm()].
#' @param ... Additional arguments passed through to [nlm()].
#'
#' @return A list containing `minimum` (minimized negative log-likelihood),
#'   `estimate` (fitted working parameter vector), `code`, `iterations`,
#'   `gradient` (`nlm()` diagnostics), `fittedModel` (fitted constrained
#'   HSMM model), `fittedHmm` (expanded HMM for the fitted model), and
#'   `workingStart` (the starting working vector).
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 200, seed = 1)
#' \donttest{
#' fit <- fitHsmm(model = testModel, dat = simData$x, stationary = TRUE)
#' fit$minimum
#' fit$fittedModel$tpm
#' }
#'
#' @export
fitHsmm <- function(
    model,
    dat,
    stationary = TRUE,
    initialType = NULL,
    workingStart = NULL,
    safeLikelihood = TRUE,
    iterlim = 500,
    stepmax = 5,
    ...
) {
  # The initial distribution used to report the fitted HMM must match the one the
  # likelihood was evaluated under, otherwise fittedHmm describes a different model
  # from the one that was optimized.
  if (is.null(initialType)) {
    initialType <- if (stationary) "stationaryHmm" else "aggregateStart"
  }
  
  if (isTRUE(model$fitShift)) {
    warning(
      "fitShift = TRUE estimates an integer shift with a gradient-based optimizer. ",
      "The likelihood is piecewise constant in the shift, so nlm() cannot move it. ",
      "Use fitHsmmShiftGrid() instead."
    )
  }
  
  if (is.null(workingStart)) {
    workingStart <- hsmmConstrainedToWorking(
      model = model,
      stationary = stationary
    )
  }
  
  objectiveFunction <- if (safeLikelihood) {
    safeHsmmNegativeLogLik
  } else {
    hsmmNegativeLogLik
  }
  
  fit <- nlm(
    f = objectiveFunction,
    p = workingStart,
    x = dat,
    modelTemplate = model,
    stationary = stationary,
    iterlim = iterlim,
    stepmax = stepmax,
    ...
  )
  
  # nlm() codes 1 and 2 indicate convergence; 3, 4 and 5 do not.
  if (fit$code >= 3) {
    warning(
      "nlm() returned convergence code ", fit$code,
      ". The fit may not be at a local optimum. ",
      "Consider raising iterlim, adjusting stepmax, or using a different start."
    )
  }

  fittedModel <- hsmmWorkingToConstrained(
    parVec = fit$estimate,
    modelTemplate = model,
    stationary = stationary
  )
  
  fittedHmm <- hsmmToHmm(
    model = fittedModel,
    initialType = initialType
  )
  
  list(
    minimum = fit$minimum,
    estimate = fit$estimate,
    code = fit$code,
    iterations = fit$iterations,
    gradient = fit$gradient,
    model = fittedModel,
    fittedModel = fittedModel,
    fittedHmm = fittedHmm,
    workingStart = workingStart,
    stationary = stationary,
    initialType = initialType
  )
}


#' Fit an HSMM by the EM algorithm
#'
#' Fits the same model as [fitHsmm()] by expectation-maximization on the
#' expanded HMM instead of direct numerical optimization. The E-step is a
#' sparse scaled forward-backward pass costing `O(nObs * (N + K^2))` (the
#' dense expanded transition matrix is never built). The M-step uses the
#' separation of the expected complete-data log-likelihood into independent
#' blocks:
#' \describe{
#'   \item{embedded TPM}{closed form, `N(i, j) / sum_j N(i, j)`}
#'   \item{`"np"`, `"npg"`, `"geo"` dwell}{closed-form life-table hazards
#'     `d_i(r) / eta_i(r)`}
#'   \item{`"pois"`, `"shiftp"`, `"nb"`, `"shiftnb"`, `"weibull"` dwell}{a few
#'     step-halved Fisher-scoring steps on the working scale (generalized EM)}
#'   \item{Normal / Poisson emissions}{closed-form weighted moments}
#'   \item{initial distribution}{closed form when `stationary = FALSE`}
#' }
#'
#' With `stationary = FALSE` this is exact (generalized) EM: the
#' log-likelihood never decreases, and fixed points are stationary points of
#' the same likelihood [fitHsmm()] maximizes. With `stationary = TRUE` the
#' initial distribution depends on every other parameter, so it is treated as
#' a plug-in and recomputed after each M-step; the result is then within
#' `O(1)` log-likelihood units of the stationary MLE rather than at it. Use
#' `polish = TRUE` to finish with a short [fitHsmm()] run from the EM
#' solution when the exact stationary MLE is needed.
#'
#' Shift parameters are held fixed, as in [fitHsmm()]; profile over them
#' with `fitHsmmShiftGrid(..., fitFunction = fitHsmmEM)`.
#'
#' Any emission family in emissionFamilies.R can be used, including the
#' functional (topological HMM) families `"bmwd"`, `"ou"` and `"nonpar"`, for
#' which `dat` is a matrix or array of curves. Decode with [hsmmViterbiEM()].
#'
#' @inheritParams fitHsmm
#' @param maxIter Maximum number of EM iterations.
#' @param tol Convergence tolerance on the relative log-likelihood change.
#' @param dwellSteps Maximum Fisher-scoring steps per EM iteration for each
#'   parametric dwell-time state.
#' @param minSd Lower bound on Normal emission standard deviations.
#' @param minCount Expected-count threshold below which a block keeps its
#'   previous value instead of being re-estimated from (near-)empty counts.
#' @param polish Logical; if `TRUE`, run [fitHsmm()] from the EM solution
#'   and return that fit, with the EM fit attached as `emFit`.
#' @param verbose Logical; print the log-likelihood each iteration.
#'
#' @return A list with the same core fields as [fitHsmm()] (`minimum`,
#'   `fittedModel`, `model`, `fittedHmm`, `stationary`, `initialType`), plus
#'   `logLik`, `iterations`, `converged`, `logLikTrace`, and `method = "EM"`.
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 300, seed = 1)
#' \donttest{
#' emFit <- fitHsmmEM(model = testModel, dat = simData$x, stationary = FALSE)
#' emFit$minimum
#' }
#'
#' @export
fitHsmmEM <- function(
    model,
    dat,
    stationary = TRUE,
    initialType = NULL,
    maxIter = 2000,
    tol = 1e-10,
    dwellSteps = 5,
    minSd = 1e-6,
    minCount = 1e-8,
    polish = FALSE,
    verbose = FALSE
) {
  if (is.null(initialType)) {
    initialType <- if (stationary) "stationaryHmm" else "aggregateStart"
  }

  if (isTRUE(model$fitShift)) {
    warning(
      "fitShift = TRUE is ignored by fitHsmmEM(): shifts are held fixed. ",
      "Use fitHsmmShiftGrid(..., fitFunction = fitHsmmEM) to profile over them."
    )
    model$fitShift <- FALSE
  }

  if (!stationary && is.null(model$delta)) {
    stop("model$delta is required when stationary = FALSE.")
  }

  control <- list(dwellSteps = dwellSteps, minSd = minSd, minCount = minCount)
  model$stationary <- stationary

  family <- getEmissionFamily(model$emissionType)

  if (polish && family$isFunctional) {
    stop("polish = TRUE needs scalar emissions: fitHsmm() has no working-scale ",
         "parametrization for functional emission parameters.")
  }

  # Reduce the data once (for curve families: per-curve sufficient statistics).
  prepared <- family$prepare(dat, model$emissionControl)
  exactEm <- !stationary && family$exactMStep(model$emissionControl)

  eStep <- hsmmEStep(model, dat, stationary, prepared)

  if (!is.finite(eStep$logLik)) {
    stop("The starting model has zero likelihood for these data.")
  }

  logLikTrace <- eStep$logLik
  converged <- FALSE
  decreases <- 0L

  for (iteration in seq_len(maxIter)) {
    candidateModel <- hsmmMStep(model, eStep, prepared, stationary, control)
    candidateE <- hsmmEStep(candidateModel, dat, stationary, prepared)

    if (!is.finite(candidateE$logLik)) {
      warning("EM step produced a zero-likelihood model; stopping at iteration ", iteration, ".")
      break
    }

    change <- candidateE$logLik - eStep$logLik

    if (change < -1e-8 * (1 + abs(eStep$logLik))) {
      decreases <- decreases + 1L

      if (exactEm) {
        warning(
          "Log-likelihood decreased by ", signif(-change, 4), " at iteration ", iteration,
          ". Exact EM is monotone, so this indicates a numerical problem."
        )
      }
    }

    model <- candidateModel
    eStep <- candidateE
    logLikTrace <- c(logLikTrace, eStep$logLik)

    if (verbose) {
      message(sprintf("EM iteration %4d: logLik = %.10f (change %.3e)", iteration, eStep$logLik, change))
    }

    if (abs(change) < tol * (1 + abs(eStep$logLik))) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning("fitHsmmEM() reached maxIter = ", maxIter, " without converging.")
  }

  emFit <- list(
    minimum = -eStep$logLik,
    logLik = eStep$logLik,
    model = model,
    fittedModel = model,
    fittedHmm = hsmmToHmm(model, initialType = initialType),
    iterations = length(logLikTrace) - 1L,
    converged = converged,
    logLikTrace = logLikTrace,
    likelihoodDecreases = decreases,
    stationary = stationary,
    initialType = initialType,
    method = "EM"
  )

  if (!polish) {
    return(emFit)
  }

  polishedFit <- fitHsmm(
    model = model,
    dat = dat,
    stationary = stationary,
    initialType = initialType
  )
  polishedFit$method <- "EM + nlm polish"
  polishedFit$emFit <- emFit
  polishedFit
}


# Index of the shift entry in a sojourn parameter vector, or NA if the family has none.
getShiftIndex <- function(sojournType) {
  if (sojournType == "shiftp") {
    return(2L)
  }
  
  if (sojournType == "shiftnb") {
    return(3L)
  }
  
  NA_integer_
}


#' Fit an HSMM by profiling over integer dwell-time shifts
#'
#' Fits a model whose dwell-time distributions include a shift parameter by
#' grid search rather than gradient-based optimization. [fitHsmm()] cannot
#' estimate a shift directly (see its documentation for why); this function
#' holds each candidate combination of shifts fixed, fits the remaining
#' parameters with [fitHsmm()] (`fitShift = FALSE`), and keeps the
#' combination with the lowest negative log-likelihood.
#'
#' @param model A constrained HSMM model, as returned by [buildHsmm()].
#' @param dat Numeric vector of observations.
#' @param shiftValues A list with one entry per HSMM state. Entries for
#'   states whose sojourn family has no shift (i.e. not `"shiftp"` or
#'   `"shiftnb"`) are ignored. A `NULL` entry keeps the model's current
#'   shift value for that state; otherwise supply a numeric vector of
#'   candidate integer shifts to try.
#' @param stationary Logical; passed through to [fitHsmm()].
#' @param initialType Passed through to [fitHsmm()].
#' @param safeLikelihood Passed through to [fitHsmm()].
#' @param maxCombinations Integer guard against an accidentally huge grid;
#'   the function errors if the shift grid exceeds this many combinations.
#' @param verbose Logical; if `TRUE`, prints progress for each combination.
#' @param fitFunction The fitter run for each shift combination:
#'   [fitHsmm()] (default) or [fitHsmmEM()].
#' @param ... Further arguments passed to `fitFunction`.
#'
#' @return The `fitHsmm()` result for the best-fitting shift combination,
#'   with two additional fields: `bestShift` (the winning shift for each
#'   state) and `shiftGrid` (a data frame of every combination tried and its
#'   negative log-likelihood).
#'
#' @examples
#' testModel <- exampleHsmmModel()
#' simData <- simulateHsmm(testModel, nObs = 200, seed = 1)
#' \donttest{
#' fit <- fitHsmmShiftGrid(
#'   model = testModel,
#'   dat = simData$x,
#'   shiftValues = list(1:2, NULL),
#'   stationary = TRUE
#' )
#' fit$bestShift
#' }
#'
#' @export
fitHsmmShiftGrid <- function(
    model,
    dat,
    shiftValues,
    stationary = TRUE,
    initialType = NULL,
    safeLikelihood = TRUE,
    maxCombinations = 64,
    verbose = FALSE,
    fitFunction = fitHsmm,
    ...
) {
  nStates <- model$nStates
  
  if (!is.list(shiftValues) || length(shiftValues) != nStates) {
    stop("shiftValues must be a list with one entry per HSMM state.")
  }
  
  shiftIndices <- vapply(model$sojournType, getShiftIndex, integer(1))
  
  candidateList <- vector("list", nStates)
  
  for (stateIndex in seq_len(nStates)) {
    shiftIndex <- shiftIndices[stateIndex]
    
    if (is.na(shiftIndex) || is.null(shiftValues[[stateIndex]])) {
      currentShift <- if (is.na(shiftIndex)) {
        NA_real_
      } else {
        model$sojournParameters[[stateIndex]][shiftIndex]
      }
      
      candidateList[[stateIndex]] <- currentShift
    } else {
      candidateList[[stateIndex]] <- as.numeric(shiftValues[[stateIndex]])
    }
  }
  
  shiftGrid <- expand.grid(candidateList, KEEP.OUT.ATTRS = FALSE)
  nCombinations <- nrow(shiftGrid)
  
  if (nCombinations > maxCombinations) {
    stop(
      "The shift grid has ", nCombinations, " combinations, above maxCombinations = ",
      maxCombinations, ". Narrow the grid or raise maxCombinations."
    )
  }
  
  fits <- vector("list", nCombinations)
  minima <- rep(NA_real_, nCombinations)
  
  # fitShift is forced off: the shift is fixed by the grid, not estimated.
  gridModelTemplate <- model
  gridModelTemplate$fitShift <- FALSE
  
  for (combinationIndex in seq_len(nCombinations)) {
    candidateModel <- gridModelTemplate
    
    for (stateIndex in seq_len(nStates)) {
      shiftIndex <- shiftIndices[stateIndex]
      
      if (!is.na(shiftIndex)) {
        candidateModel$sojournParameters[[stateIndex]][shiftIndex] <-
          shiftGrid[combinationIndex, stateIndex]
      }
    }
    
    fitArguments <- list(
      model = candidateModel,
      dat = dat,
      stationary = stationary,
      initialType = initialType,
      ...
    )

    # safeLikelihood only means something to the nlm()-based fitter.
    if ("safeLikelihood" %in% names(formals(fitFunction))) {
      fitArguments$safeLikelihood <- safeLikelihood
    }

    candidateFit <- tryCatch(
      do.call(fitFunction, fitArguments),
      error = function(errorCondition) NULL
    )
    
    fits[[combinationIndex]] <- candidateFit
    
    if (!is.null(candidateFit)) {
      minima[combinationIndex] <- candidateFit$minimum
    }
    
    if (verbose) {
      message(
        "shift combination ", combinationIndex, " of ", nCombinations,
        ": nll = ", signif(minima[combinationIndex], 8)
      )
    }
  }
  
  if (all(is.na(minima))) {
    stop("Every shift combination failed to fit.")
  }
  
  bestIndex <- which.min(minima)
  
  gridSummary <- cbind(
    shiftGrid,
    negativeLogLik = minima
  )
  
  names(gridSummary)[seq_len(nStates)] <- paste0("shiftState", seq_len(nStates))
  
  output <- fits[[bestIndex]]
  output$shiftGrid <- gridSummary
  output$bestShift <- as.numeric(shiftGrid[bestIndex, ])
  
  output
}
