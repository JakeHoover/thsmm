# Refined simulation studies for the THSMM paper, with the CRAN THMM package.
#
# Run from the package root:
#   Rscript tests/paperSimsRefined.R            all sections
#   Rscript tests/paperSimsRefined.R 1 3        selected sections
#   Rscript tests/paperSimsRefined.R --quick    2 replicates, fewer curves: a smoke test
# The number of replicates and the other run settings are set in the "Run
# settings" block below: start with nReps = 1 for a proof of concept, then
# switch to 50. Environment variables override that block when set:
#   PAPER_SIMS_REPS    replicates per setting
#   PAPER_SIMS_CORES   parallel workers
#
# Compared with tests/paperSims.R, the scenarios are fewer and chosen to cover
# distinct situations: the Section 3 showcases are the four settings in which a
# pilot run (8 replicates each, 24 candidates) found the largest THSMM gains
# over both THMMs (bimodal, long minimum, near-deterministic and heterogeneous
# dwell times; the peaked negative binomial of paperSims.R gained little and
# was dropped). Sensitivity blocks vary emission separation and sequence
# length, and Section 2 adds a sample-size series.
#
# Emission settings (Kashlak, Loliencar and Heo, 2023, Sec. 5), curves on 101
# equispaced points of [0, 1] with O(0) = 0:
#   BM    Brownian motion with linear drift c           ("bmwd", hurst 1/2)
#   OU    Ornstein-Uhlenbeck dO = (b0 - b1 O) dtau + dW  ("ou", b0 = c mu, b1 = c)
#   fBM   fractional BM, hurst 0.75, linear drift c      ("bmwd", hurst known)
#
# THMM fits (geometric dwell times), all from the same starting point:
#   THMM (funHMM)       the CRAN package funHMM (Kashlak), thmm() with its
#                       automatic sigma. It follows the paper's Onsager-Machlup
#                       functionals but discretizes them differently from this
#                       package (BM: V = 1 rather than the observed span; OU:
#                       Stratonovich form), so it is close to, not identical
#                       with, the package's geometric fit
#   THSMM, geometric    the THMM as a special case of this package
# funHMM gets the curves without the O(0) column: 100 points at k / 100, its
# own data convention (as funHMM::rthmm() simulates them).
#
# THSMM fits:
#   THSMM, npg(12)      nonparametric dwell head (11 free probabilities) +
#                       geometric tail: no knowledge of the dwell family
#   THSMM, true family  the generating family per state; integer shifts are
#                       profiled over a grid by fitHsmmShiftGrid()
#
# Section 1 -- THMM recovery (geometric truth): baseline, low separation,
# sticky chain (long dwells), fast-switching chain, and the paper's
# five-state A1 setting. Does the THSMM reproduce the THMM fits, and does the
# npg THSMM find the geometric shape on its own?
#
# Section 2 -- sojourn recovery (medium separation, 1000 curves): long minimum
# dwell (shifted Poisson), overdispersed (negative binomial), near-
# deterministic (Weibull) and bimodal (nonparametric) dwell times, plus a
# sample-size series (250 to 2000 curves, negative binomial, BM). The THMM
# fits show what a geometric dwell implies for the same data.
#
# Section 3 -- Viterbi ARI, THMM vs THSMM, non-geometric truth at low
# separation (500 curves): bimodal, long minimum, near-deterministic and
# heterogeneous dwell times, plus sensitivity series for the bimodal truth
# (BM) over emission separation (very low / low / medium) and sequence length
# (200 / 500 / 1000 curves). Paired differences per data set.
#
# Output in tests/output/paperSimsRefined/ (tests/output/paperSimsRefined-poc/
# when nReps = 1): per-replicate results
# (sectionN.rds, sectionN-replicates.csv), summary tables (.csv and LaTeX
# .tex) and figures (sectionN.pdf).

# ---------------------------------------------------------------------------
# Run settings: edit these
# ---------------------------------------------------------------------------

nReps <- 1          # replicates per setting: 1 = proof of concept, 50 for the paper
nCores <- parallel::detectCores() - 1

# A proof-of-concept run (nReps = 1) writes to its own folder, so it never
# overwrites the results of a full run.
outputName <- if (nReps == 1) "paperSimsRefined-poc" else "paperSimsRefined"

# ---------------------------------------------------------------------------

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

if (!requireNamespace("funHMM", quietly = TRUE)) {
  stop("funHMM is not installed: install.packages(\"funHMM\")")
}

library(parallel)

arguments <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% arguments
sections <- as.integer(setdiff(arguments, "--quick"))
if (length(sections) == 0L) sections <- 1:3

nReps <- if (quick) 2L else as.integer(Sys.getenv("PAPER_SIMS_REPS", nReps))
nCores <- max(1L, as.integer(Sys.getenv("PAPER_SIMS_CORES", nCores)))
if (Sys.getenv("PAPER_SIMS_REPS") != "") outputName <- if (nReps == 1) "paperSimsRefined-poc" else "paperSimsRefined"

outputDir <- file.path("tests", "output", outputName)
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

grid <- seq(0, 1, length.out = 101)
curveScale <- if (quick) 0.4 else 1   # quick mode shortens every sequence
npgSize <- 12                     # npg: 11 free head probabilities + geometric tail
parametricDwellSize <- 30         # exact dwell support before the geometric tail
maxDwell <- 300                   # support for mean dwells and TV distances
maxPlot <- 20                     # dwell support shown in the figures
emTol <- 1e-8
emMaxIter <- 2000

ouParams <- function(mu, rate) cbind(b0 = rate * mu, b1 = rate)

# Emission parameters by separation level. Single-curve classification at the
# true parameters is roughly 0.70 ("vlow"), 0.78 ("low") and 0.9+ ("medium")
# for three states; "paper5" are five-state settings of Kashlak et al. (2023,
# Sec. 5). States are ordered by drift / mean level, as initializeHsmmEM()
# orders its starting states.
processes <- list(
  BM = list(
    label = "Brownian motion with drift", emissionType = "bmwd", hurst = 0.5, control = list(grid = grid),
    params = list(vlow = cbind(drift = c(-1.5, 0, 1.5)), low = cbind(drift = c(-2, 0, 2)),
                  medium = cbind(drift = c(-3, 0, 3)), paper5 = cbind(drift = c(-4, -2, 0, 2, 4)))
  ),
  OU = list(
    label = "Ornstein-Uhlenbeck", emissionType = "ou", hurst = 0.5, control = list(grid = grid),
    params = list(vlow = ouParams(c(-1, 0, 1), 1.5), low = ouParams(c(-1, 0, 1), 2), medium = ouParams(c(-1, 0, 1), 3),
                  paper5 = ouParams(c(-2, 0, 1, 2, 4), c(4, 4, 20, 2, 8)))   # the paper's (mu, c), ordered by mu
  ),
  fBM = list(
    label = "Fractional BM with drift (H = 0.75)", emissionType = "bmwd", hurst = 0.75,
    control = list(grid = grid, hurst = 0.75),
    params = list(vlow = cbind(drift = c(-1.5, 0, 1.5)), low = cbind(drift = c(-2, 0, 2)),
                  medium = cbind(drift = c(-3, 0, 3)), paper5 = cbind(drift = c(-4, -2, 0, 2, 4)))
  )
)
allProcesses <- names(processes)

withSize <- function(...) c(..., dwellSize = parametricDwellSize)
geometric <- function(leave) list(type = rep("geo", length(leave)), params = lapply(leave, function(p) c(prob = p)))
bimodal <- function(first, second, weight = 0.5, size = 12) {
  pmf <- numeric(size)
  pmf[first] <- weight
  pmf[second] <- 1 - weight
  pmf
}

# A scenario: dwell types and parameters per state, separation level,
# number of curves, the processes it runs on, and for sensitiI'm excited to see the results. vity series a
# `series` name and the `value` plotted on its x-axis.
scenario <- function(dwell, level, nCurves, processes = allProcesses, shiftCandidates = NULL, series = NA, value = NA) {
  list(type = dwell$type, params = dwell$params, level = level, nCurves = round(nCurves * curveScale),
       processes = processes, shiftCandidates = shiftCandidates, series = series, value = value)
}

# Section 1: geometric truth.
section1Scenarios <- list(
  Baseline = scenario(geometric(c(0.4, 0.25, 0.15)), "medium", 500),          # mean dwells 2.5, 4, 6.7
  `Low separation` = scenario(geometric(c(0.4, 0.25, 0.15)), "low", 500),
  `Sticky chain` = scenario(geometric(c(0.05, 0.08, 0.12)), "medium", 1000),  # mean dwells 20, 12.5, 8.3
  `Fast switching` = scenario(geometric(c(0.8, 0.7, 0.6)), "medium", 500),    # mean dwells 1.25, 1.4, 1.7
  # The paper's A1: stay 0.64, move 0.09 to each other state.
  `Five states (paper A1)` = scenario(geometric(rep(0.36, 5)), "paper5", 500)
)

# Section 2: non-geometric dwell families, medium separation.
longMinimumDwell <- list(type = rep("shiftp", 3),
                         params = list(withSize(lambda = 1, shift = 6), withSize(lambda = 2, shift = 7),
                                       withSize(lambda = 1, shift = 8)))                     # means 7, 9, 9
overdispersed <- list(type = rep("nb", 3),
                      params = list(withSize(size = 2, prob = 0.4), withSize(size = 5, prob = 0.5),
                                    withSize(size = 10, prob = 0.55)))                       # 1 + NB: means 4, 6, 9.2
nearDeterministic <- list(type = rep("weibull", 3),
                          params = list(withSize(shape = 6, scale = 5), withSize(shape = 6, scale = 6),
                                        withSize(shape = 6, scale = 7)))                     # means 4.9, 6.0, 7.0
bimodalDwell <- list(type = rep("np", 3), params = list(bimodal(2, 8), bimodal(3, 9), bimodal(2, 10)))
heterogeneous <- list(type = c("pois", "shiftp", "nb"),                                       # short, long, peaked
                      params = list(withSize(lambda = 1), withSize(lambda = 1, shift = 8),
                                    withSize(size = 20, prob = 20 / 24)))                    # means 2, 9, 5

section2Scenarios <- list(
  `Long minimum dwell` = scenario(longMinimumDwell, "medium", 1000, shiftCandidates = 6:9),
  Overdispersed = scenario(overdispersed, "medium", 1000),
  `Near-deterministic` = scenario(nearDeterministic, "medium", 1000),
  Bimodal = scenario(bimodalDwell, "medium", 1000)
)
for (n in c(250, 500, 2000)) {
  section2Scenarios[[sprintf("Overdispersed, %d curves", n)]] <-
    scenario(overdispersed, "medium", n, processes = "BM", series = "Sample size", value = n)
}

# Section 3: the showcases, low separation, 500 curves.
section3Scenarios <- list(
  Bimodal = scenario(bimodalDwell, "low", 500),
  `Long minimum dwell` = scenario(longMinimumDwell, "low", 500, shiftCandidates = 6:9),
  `Near-deterministic` = scenario(nearDeterministic, "low", 500),
  # The true family of each state is assigned in drift / mean-level order (the
  # order of both the truth and initializeHsmmEM()'s starting states).
  Heterogeneous = scenario(heterogeneous, "low", 500, shiftCandidates = 6:9)
)
sensitivity <- list(
  `Bimodal, very low separation` = scenario(bimodalDwell, "vlow", 500, "BM", series = "Separation", value = 1),
  `Bimodal, medium separation` = scenario(bimodalDwell, "medium", 500, "BM", series = "Separation", value = 3),
  `Bimodal, 200 curves` = scenario(bimodalDwell, "low", 200, "BM", series = "Sequence length", value = 200),
  `Bimodal, 1000 curves` = scenario(bimodalDwell, "low", 1000, "BM", series = "Sequence length", value = 1000)
)
section3Scenarios <- c(section3Scenarios, sensitivity)

if (quick) {
  # Only the true shifts: the grid is the costliest part of a fit.
  for (s in c("section1Scenarios", "section2Scenarios", "section3Scenarios")) {
    scenarios <- get(s)
    for (name in names(scenarios)) {
      if (!is.null(scenarios[[name]]$shiftCandidates)) {
        shifts <- unlist(lapply(scenarios[[name]]$params, function(p) if ("shift" %in% names(p)) p[["shift"]]))
        scenarios[[name]]$shiftCandidates <- sort(unique(shifts))
      }
    }
    assign(s, scenarios)
  }
}


# ---------------------------------------------------------------------------
# Helpers: models and fits
# ---------------------------------------------------------------------------

permutationCache <- list()
permutations <- function(n) {
  key <- as.character(n)
  if (!is.null(permutationCache[[key]])) return(permutationCache[[key]])
  build <- function(n) {
    if (n == 1L) return(matrix(1L, 1, 1))
    smaller <- build(n - 1L)
    do.call(rbind, lapply(seq_len(n), function(first) cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))))
  }
  permutationCache[[key]] <<- build(n)
}

# match[i] = estimated state that best matches true state i (maximizes matched counts).
matchStates <- function(estimated, truth, nStates) {
  tab <- table(factor(estimated, seq_len(nStates)), factor(truth, seq_len(nStates)))
  all <- permutations(nStates)
  all[which.max(apply(all, 1, function(p) sum(tab[cbind(p, seq_len(nStates))]))), ]
}

adjustedRand <- function(a, b) funHMM::ari(a, b)

seedFor <- function(section, processIndex, scenarioIndex, replicate) {
  as.integer(1e6 * section + 1e4 * processIndex + 1e2 * scenarioIndex + replicate)
}

makeTruth <- function(process, scenario) {
  nStates <- length(scenario$type)
  omega <- matrix(1 / (nStates - 1), nStates, nStates)
  diag(omega) <- 0
  buildHsmm(omega, scenario$params, scenario$type, process$params[[scenario$level]], process$emissionType,
            delta = rep(1 / nStates, nStates), stationary = FALSE, emissionControl = process$control)
}

dwellSizeFor <- function(sojournType, scenario) {
  ifelse(sojournType == "npg", npgSize,
         ifelse(sojournType == "np", vapply(scenario$params, length, integer(1))[seq_along(sojournType)],
                parametricDwellSize))
}

# Starting model from the data alone; deterministic given the seed, so every
# dwell family shares the same emission start.
startModel <- function(x, process, sojournType, seed, scenario) {
  nStates <- length(scenario$type)
  sojournType <- rep_len(sojournType, nStates)
  initializeHsmmEM(x, nStates, sojournType, process$emissionType, process$control,
                   dwellSize = dwellSizeFor(sojournType, scenario), stationary = FALSE, seed = seed)$model
}

fitThsmm <- function(start, x, scenario) {
  if (any(start$sojournType %in% c("shiftp", "shiftnb"))) {
    # Integer shifts are held fixed by EM: profile them over the scenario's grid.
    shiftValues <- lapply(start$sojournType, function(type) if (type %in% c("shiftp", "shiftnb")) scenario$shiftCandidates)
    return(suppressWarnings(fitHsmmShiftGrid(start, x, shiftValues, stationary = FALSE, fitFunction = fitHsmmEM,
                                             tol = emTol, maxIter = emMaxIter, maxCombinations = 256)))
  }
  suppressWarnings(fitHsmmEM(start, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
}

# THMM transition matrix of a geometric-dwell model: A = diag(1 - p) + diag(p) omega.
thmmMatrix <- function(model) {
  leave <- vapply(model$sojournParameters, function(s) s[["prob"]], numeric(1))
  diag(1 - leave) + leave * model$tpm
}

# Package emission parameters <-> funHMM's (drift vector, or (mean, rate) for OU).
toFunHmmParams <- function(params, emissionType) {
  if (emissionType == "ou") {
    rate <- pmax(params[, 2], 1e-3)   # a start at the rate = 0 boundary has no mean
    return(cbind(mean = params[, 1] / rate, rate = rate))
  }
  params[, 1]
}
fromFunHmmParams <- function(par, emissionType) {
  if (emissionType == "ou") {
    par <- matrix(par, ncol = 2)
    return(cbind(b0 = par[, 1] * par[, 2], b1 = par[, 2]))   # b0 is NA when the rate is 0
  }
  cbind(drift = as.vector(par))
}

fitFunHmm <- function(x, start, process) {
  fit <- suppressWarnings(funHMM::thmm(
    x[, -1], start$nStates, type = process$emissionType, hurst = process$hurst,
    par = toFunHmmParams(start$emissionParameters, process$emissionType),
    trans = thmmMatrix(start), init.prob = start$delta, tol = emTol, max.iter = emMaxIter
  ))
  list(states = fit$states, A = fit$trans, emissionParameters = fromFunHmmParams(fit$par, process$emissionType),
       extra = c(sigma = fit$sigma, converged = as.numeric(fit$converged)))
}


# ---------------------------------------------------------------------------
# Helpers: scoring
# ---------------------------------------------------------------------------

# Fitted methods in a common form: decoded states, emission parameters, and
# per-(estimated)-state dwell pmf / family / parameters.
fromPackageFit <- function(fit, x) {
  model <- fit$fittedModel
  list(
    states = hsmmViterbiEM(model, x)$states,
    emissionParameters = model$emissionParameters,
    dwellPmf = function(j) getModelDwellPmf(model, j, maxDwell),
    dwellType = function(j) model$sojournType[j],
    dwellParams = function(j) model$sojournParameters[[j]],
    A = if (all(model$sojournType == "geo")) thmmMatrix(model)
  )
}

fromThmmFit <- function(fit) {
  if (is.null(fit)) return(NULL)
  c(fit, list(
    dwellPmf = function(j) dgeom(seq_len(maxDwell) - 1, 1 - fit$A[j, j]),
    dwellType = function(j) "geo",
    dwellParams = function(j) c(prob = 1 - fit$A[j, j])
  ))
}

emissionSummary <- function(params, emissionType) {
  if (emissionType == "ou") return(cbind(b0 = params[, 1], b1 = params[, 2], mu = params[, 1] / params[, 2]))
  cbind(drift = params[, 1])
}

resultRow <- function(state, parameter, truth, estimate) {
  data.frame(state = state, parameter = parameter, truth = truth, estimate = estimate)
}

# Score one method: ARI, emission and dwell parameters of matched states, mean
# dwell, TV distance of the dwell pmf, THMM transition-matrix error (geometric
# truth), plus the matched dwell pmfs on 1:maxPlot for figures.
scoreMethod <- function(method, fitted, truth, trueStates) {
  nStates <- truth$nStates
  match <- matchStates(fitted$states, trueStates, nStates)
  rows <- list(resultRow(NA, "ARI", 1, adjustedRand(fitted$states, trueStates)))
  if (!is.null(fitted$extra)) {
    for (p in names(fitted$extra)) rows[[length(rows) + 1L]] <- resultRow(NA, p, NA, fitted$extra[[p]])
  }

  trueEmission <- emissionSummary(truth$emissionParameters, truth$emissionType)
  estEmission <- emissionSummary(fitted$emissionParameters, truth$emissionType)[match, , drop = FALSE]
  dwellTimes <- seq_len(maxDwell)
  pmfRows <- list()

  for (i in seq_len(nStates)) {
    j <- match[i]
    for (p in colnames(trueEmission)) rows[[length(rows) + 1L]] <- resultRow(i, p, trueEmission[i, p], estEmission[i, p])

    truePmf <- getModelDwellPmf(truth, i, maxDwell)
    estPmf <- fitted$dwellPmf(j)
    rows[[length(rows) + 1L]] <- resultRow(i, "mean dwell", sum(dwellTimes * truePmf), sum(dwellTimes * estPmf))
    rows[[length(rows) + 1L]] <- resultRow(i, "TV dwell", 0, 0.5 * sum(abs(truePmf - estPmf)))

    # Dwell parameters, when the fitted family is the true one (dwellSize is fixed).
    if (fitted$dwellType(j) == truth$sojournType[i] && truth$sojournType[i] != "np") {
      trueDwell <- truth$sojournParameters[[i]]
      estDwell <- fitted$dwellParams(j)
      for (p in setdiff(names(trueDwell), "dwellSize")) {
        rows[[length(rows) + 1L]] <- resultRow(i, p, trueDwell[[p]], estDwell[[p]])
      }
    }

    pmfRows[[i]] <- data.frame(method = method, state = i, dwell = seq_len(maxPlot), prob = estPmf[seq_len(maxPlot)])
  }

  if (!is.null(fitted$A) && all(truth$sojournType == "geo")) {
    rows[[length(rows) + 1L]] <- resultRow(NA, "max |A error|", 0, max(abs(fitted$A[match, match] - thmmMatrix(truth))))
  }

  list(rows = cbind(method = method, do.call(rbind, rows)), pmf = do.call(rbind, pmfRows))
}

# Run the methods (named list of functions returning the common form, or NULL
# on failure), score each, and record the agreement of every other method's
# Viterbi path with the funHMM THMM.
scoreMethods <- function(methods, truth, sim) {
  fitted <- list()
  scores <- list()
  for (method in names(methods)) {
    elapsed <- system.time(result <- tryCatch(methods[[method]](), error = function(e) NULL))[["elapsed"]]
    if (is.null(result)) {
      scores[[method]] <- list(rows = data.frame(method = method, state = NA, parameter = "failed", truth = NA, estimate = 1))
      next
    }
    fitted[[method]] <- result
    scores[[method]] <- scoreMethod(method, result, truth, sim$states)
    scores[[method]]$rows <- rbind(scores[[method]]$rows,
                                   data.frame(method = method, state = NA, parameter = "seconds", truth = NA, estimate = elapsed))
  }

  reference <- fitted[["THMM (funHMM)"]]
  agreement <- if (!is.null(reference)) do.call(rbind, lapply(setdiff(names(fitted), "THMM (funHMM)"), function(m) {
    data.frame(method = m, state = NA, parameter = "ARI vs funHMM", truth = NA,
               estimate = adjustedRand(fitted[[m]]$states, reference$states))
  }))

  list(rows = rbind(do.call(rbind, lapply(scores, `[[`, "rows")), agreement),
       pmf = do.call(rbind, lapply(scores, `[[`, "pmf")))
}

# The methods of one replicate. `which` picks THMM and THSMM fits.
methodList <- function(x, process, scenario, seed, which) {
  geoStart <- startModel(x, process, "geo", seed, scenario)
  methods <- list(
    "THMM (funHMM)" = function() fromThmmFit(fitFunHmm(x, geoStart, process)),
    "THSMM, geometric" = function() fromPackageFit(fitThsmm(geoStart, x, scenario), x),
    "THSMM, npg(12)" = function() fromPackageFit(fitThsmm(startModel(x, process, "npg", seed, scenario), x, scenario), x),
    "THSMM, true family" = function() {
      fromPackageFit(fitThsmm(startModel(x, process, scenario$type, seed, scenario), x, scenario), x)
    }
  )
  methods[which]
}

replicateFunction <- function(section, scenarios, methods) {
  function(task) {
    process <- processes[[task$process]]
    scenario <- scenarios[[task$scenario]]
    seed <- seedFor(section, match(task$process, names(processes)), match(task$scenario, names(scenarios)), task$replicate)
    truth <- makeTruth(process, scenario)
    sim <- simulateHsmm(truth, scenario$nCurves, initialType = "aggregateStart", seed = seed)
    scoreMethods(methodList(sim$x, process, scenario, seed, methods), truth, sim)
  }
}

makeTasks <- function(scenarios) {
  do.call(rbind, lapply(names(scenarios), function(name) {
    expand.grid(process = scenarios[[name]]$processes, scenario = name, replicate = seq_len(nReps), stringsAsFactors = FALSE)
  }))
}

runTasks <- function(tasks, replicateFunction, section) {
  started <- Sys.time()
  outputs <- mclapply(seq_len(nrow(tasks)), function(t) {
    tryCatch(replicateFunction(tasks[t, , drop = FALSE]), error = function(e) list(error = conditionMessage(e)))
  }, mc.cores = nCores, mc.preschedule = FALSE)

  failed <- vapply(outputs, function(o) !is.list(o) || !is.null(o$error), logical(1))
  if (any(failed)) {
    example <- outputs[[which(failed)[1]]]
    cat(sprintf("  %d of %d replicates failed, e.g.: %s\n", sum(failed), length(outputs),
                if (is.list(example)) example$error else as.character(example)))
  }

  attach <- function(field) do.call(rbind, lapply(which(!failed), function(t) {
    if (is.null(outputs[[t]][[field]])) return(NULL)
    cbind(tasks[t, , drop = FALSE], outputs[[t]][[field]], row.names = NULL)
  }))
  cat(sprintf("  Section %d: %d replicates in %.1f minutes on %d cores\n", section, nrow(tasks),
              as.numeric(difftime(Sys.time(), started, units = "mins")), nCores))
  list(rows = attach("rows"), pmf = attach("pmf"))
}


# ---------------------------------------------------------------------------
# Helpers: summaries and output
# ---------------------------------------------------------------------------

methodOrder <- c("THMM (funHMM)", "THSMM, geometric", "THSMM, npg(12)", "THSMM, true family")
methodColours <- c("THMM (funHMM)" = "#000000", "THSMM, geometric" = "#0072B2",
                   "THSMM, npg(12)" = "#D55E00", "THSMM, true family" = "#009E73")   # Okabe-Ito

# Bias / sd / RMSE of every (method, state, parameter) across replicates.
summarizeEstimates <- function(rows, by) {
  rows <- rows[!is.na(rows$truth) & !rows$parameter %in% c("seconds", "ARI"), ]
  groups <- split(rows, rows[, by], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    error <- g$estimate - g$truth
    cbind(g[1, by, drop = FALSE],
          data.frame(truth = g$truth[1], mean = mean(g$estimate, na.rm = TRUE), bias = mean(error, na.rm = TRUE),
                     sd = sd(g$estimate, na.rm = TRUE), rmse = sqrt(mean(error^2, na.rm = TRUE)),
                     n = sum(!is.na(g$estimate))))
  }))
  rownames(out) <- NULL
  out[order(out$process, out$scenario, match(out$method, methodOrder)), ]
}

meanSd <- function(v, digits = 3) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_character_)
  sprintf("%.*f (%.*f)", digits, mean(v), digits, sd(v))
}

# One row per (process, scenario, method): mean (sd) of the listed metrics.
headlineTable <- function(rows, metrics) {
  groups <- split(rows, rows[, c("process", "scenario", "method")], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    cells <- lapply(metrics, function(p) meanSd(g$estimate[g$parameter == p]))
    names(cells) <- metrics
    cbind(g[1, c("process", "scenario", "method")], as.data.frame(cells, check.names = FALSE),
          failures = sum(g$parameter == "failed"))
  }))
  rownames(out) <- NULL
  out <- out[, c(TRUE, TRUE, TRUE, vapply(out[, -(1:3), drop = FALSE], function(column) any(!is.na(column)), logical(1)))]
  out[order(match(out$process, allProcesses), out$scenario, match(out$method, methodOrder)), ]
}

# Paired differences (same data set) of each THSMM against each THMM fit.
pairedTable <- function(rows, metric, higherIsBetter = TRUE) {
  values <- rows[rows$parameter == metric, c("process", "scenario", "replicate", "method", "estimate")]
  if (metric == "TV dwell") values <- aggregate(estimate ~ process + scenario + replicate + method, values, mean)
  wide <- reshape(values, idvar = c("process", "scenario", "replicate"), timevar = "method", direction = "wide")
  names(wide) <- sub("^estimate\\.", "", names(wide))
  thsmms <- intersect(c("THSMM, npg(12)", "THSMM, true family"), names(wide))
  thmms <- intersect(c("THMM (funHMM)", "THSMM, geometric"), names(wide))
  settings <- unique(wide[, c("process", "scenario")])
  out <- do.call(rbind, lapply(seq_len(nrow(settings)), function(s) {
    g <- merge(settings[s, ], wide)
    do.call(rbind, lapply(thsmms, function(a) do.call(rbind, lapply(thmms, function(b) {
      difference <- g[[a]] - g[[b]]
      if (!higherIsBetter) difference <- -difference
      difference <- difference[!is.na(difference)]
      cbind(settings[s, ], data.frame(
        comparison = paste(a, "vs", sub("THSMM, geometric", "THMM (package)", b)),
        `mean gain` = mean(difference), se = sd(difference) / sqrt(length(difference)),
        `share THSMM better` = mean(difference > 0),
        `Wilcoxon p` = if (length(difference) > 1 && any(difference != 0)) suppressWarnings(wilcox.test(difference)$p.value) else NA,
        n = length(difference), check.names = FALSE))
    }))))
  }))
  rownames(out) <- NULL
  out[order(match(out$process, allProcesses), out$scenario), ]
}

writeLatex <- function(table, file, caption, digits = 3) {
  escape <- function(s) gsub("([&%#_$])", "\\\\\\1", s)
  cells <- vapply(table, function(column) {
    if (is.numeric(column)) formatC(column, digits = digits, format = "g") else escape(as.character(column))
  }, character(nrow(table)))
  cells <- matrix(cells, nrow = nrow(table))
  align <- paste(ifelse(vapply(table, is.numeric, logical(1)), "r", "l"), collapse = "")
  lines <- c("\\begin{table}[ht]", "\\centering", sprintf("\\caption{%s}", caption),
             sprintf("\\begin{tabular}{%s}", align), "\\hline",
             paste(paste(escape(names(table)), collapse = " & "), "\\\\"), "\\hline",
             paste(apply(cells, 1, paste, collapse = " & "), "\\\\"),
             "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
}

saveSection <- function(section, result, tables) {
  saveRDS(result, file.path(outputDir, sprintf("section%d.rds", section)))
  write.csv(result$rows, file.path(outputDir, sprintf("section%d-replicates.csv", section)), row.names = FALSE)
  for (name in names(tables)) {
    write.csv(tables[[name]]$table, file.path(outputDir, sprintf("section%d-%s.csv", section, name)), row.names = FALSE)
    writeLatex(tables[[name]]$table, file.path(outputDir, sprintf("section%d-%s.tex", section, name)), tables[[name]]$caption)
  }
}

printTable <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(format(table, digits = 3), row.names = FALSE)
}

# Dwell pmfs: truth (grey bars), fitted median per method (lines) and the
# 10-90% band of the methods listed in `bands`.
plotDwellPanel <- function(pmf, truth, i, main, bands) {
  truePmf <- getModelDwellPmf(truth, i, maxPlot)
  methods <- intersect(methodOrder, unique(pmf$method))
  quantiles <- lapply(methods, function(m) {
    d <- pmf[pmf$method == m & pmf$state == i, ]
    sapply(split(d$prob, d$dwell), quantile, c(0.1, 0.5, 0.9), na.rm = TRUE)
  })
  names(quantiles) <- methods
  upper <- max(truePmf, unlist(lapply(quantiles, function(q) q[3, ])))
  plot(seq_len(maxPlot), truePmf, type = "h", lwd = 6, col = "grey80", ylim = c(0, upper * 1.05),
       xlab = "dwell time (curves)", ylab = "probability", main = main, font.main = 1, cex.main = 0.95)
  for (m in intersect(bands, methods)) {
    polygon(c(seq_len(maxPlot), rev(seq_len(maxPlot))), c(quantiles[[m]][1, ], rev(quantiles[[m]][3, ])),
            col = adjustcolor(methodColours[[m]], 0.2), border = NA)
  }
  for (m in methods) lines(seq_len(maxPlot), quantiles[[m]][2, ], type = "b", pch = 19, cex = 0.45, col = methodColours[[m]])
}

methodLegend <- function(position, methods, truth = TRUE, ...) {
  legend(position, c(if (truth) "truth", methods), col = c(if (truth) "grey80", methodColours[methods]),
         lwd = c(if (truth) 6, rep(1.5, length(methods))), bty = "n", ...)
}

# Boxplots of one metric by method, for each scenario of a process.
metricBoxplots <- function(rows, metric, scenarios, ylab, title, ylim = NULL) {
  for (name in allProcesses) {
    shown <- names(scenarios)[vapply(scenarios, function(s) name %in% s$processes && is.na(s$series), logical(1))]
    if (length(shown) == 0L) next
    par(mfrow = c(1, length(shown)), mar = c(3, 4, 2.5, 0.5), oma = c(2.5, 0, 2, 0), las = 1)
    for (scenarioName in shown) {
      g <- rows[rows$process == name & rows$scenario == scenarioName & rows$parameter == metric, ]
      if (metric == "TV dwell") g <- aggregate(estimate ~ replicate + method, g, mean)
      methods <- intersect(methodOrder, unique(g$method))
      boxplot(estimate ~ factor(method, methods), data = g, col = adjustcolor(methodColours[methods], 0.6),
              names = seq_along(methods), xlab = "", ylab = ylab, main = scenarioName, font.main = 1,
              ylim = ylim %||% range(g$estimate, na.rm = TRUE), outcex = 0.5)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", paste(seq_along(methods), methods), fill = adjustcolor(methodColours[methods], 0.6),
           horiz = TRUE, bty = "n", cex = 0.75)
    mtext(paste(title, processes[[name]]$label, sep = ": "), side = 3, line = -1.5, font = 2)
  }
}

cat(sprintf("paperSimsRefined: sections %s, %d replicates, %d cores%s\n", paste(sections, collapse = ", "),
            nReps, nCores, if (quick) " (quick mode)" else ""))


# ===========================================================================
# Section 1: can the THSMM recover the THMM?
# ===========================================================================

if (1 %in% sections) {
  cat("\nSection 1: THMM recovery\n")
  methods1 <- c("THMM (funHMM)", "THSMM, geometric", "THSMM, npg(12)")
  result1 <- runTasks(makeTasks(section1Scenarios), replicateFunction(1, section1Scenarios, methods1), 1)
  rows <- result1$rows

  headline <- headlineTable(rows, c("ARI", "ARI vs funHMM", "TV dwell", "max |A error|", "sigma", "seconds"))
  estimates <- summarizeEstimates(rows, c("process", "scenario", "method", "state", "parameter"))
  estimates <- estimates[estimates$parameter != "TV dwell", ]

  saveSection(1, result1, list(
    summary = list(table = headline, caption = paste(
      "THMM recovery (geometric truth): Viterbi ARI against the truth and against the funHMM THMM,",
      "TV distance between fitted and true dwell pmfs, max abs error of the THMM transition matrix,",
      "and funHMM's automatic sigma; mean (sd) over replicates.")),
    estimates = list(table = estimates, caption = "THMM recovery: parameter estimates of matched states.")
  ))
  printTable("THMM recovery: mean (sd) over replicates", headline)

  pdf(file.path(outputDir, "section1.pdf"), width = 11, height = 4)
  metricBoxplots(rows, "ARI", section1Scenarios, "Viterbi ARI vs truth", "Section 1, THMM recovery")
  metricBoxplots(rows, "ARI vs funHMM", section1Scenarios, "ARI vs funHMM path", "Section 1, agreement with funHMM",
                 ylim = c(0, 1))
  # Does the npg THSMM find the geometric shape? Baseline scenario, every process.
  for (name in allProcesses) {
    truth <- makeTruth(processes[[name]], section1Scenarios$Baseline)
    p <- result1$pmf[result1$pmf$process == name & result1$pmf$scenario == "Baseline", ]
    par(mfrow = c(1, truth$nStates), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0), las = 1)
    for (i in seq_len(truth$nStates)) {
      plotDwellPanel(p, truth, i, sprintf("State %d: geometric, p = %.2f", i, truth$sojournParameters[[i]][["prob"]]),
                     "THSMM, npg(12)")
      if (i == 1) methodLegend("topright", intersect(methodOrder, unique(p$method)), cex = 0.7)
    }
    mtext(paste("Section 1, baseline dwell pmfs:", processes[[name]]$label), outer = TRUE, font = 2)
  }
  invisible(dev.off())
}


# ===========================================================================
# Section 2: recovery of sojourn distributions
# ===========================================================================

if (2 %in% sections) {
  cat("\nSection 2: sojourn recovery\n")
  methods2 <- c("THMM (funHMM)", "THSMM, geometric", "THSMM, npg(12)", "THSMM, true family")
  result2 <- runTasks(makeTasks(section2Scenarios), replicateFunction(2, section2Scenarios, methods2), 2)
  rows <- result2$rows

  headline <- headlineTable(rows, c("ARI", "TV dwell", "seconds"))
  estimates <- summarizeEstimates(rows, c("process", "scenario", "method", "state", "parameter"))
  dwellParameters <- c("lambda", "shift", "size", "prob", "shape", "scale", "mean dwell", "TV dwell")
  dwellTable <- estimates[estimates$parameter %in% dwellParameters, ]
  emissionTable <- estimates[!estimates$parameter %in% dwellParameters, ]
  pairedTv <- pairedTable(rows, "TV dwell", higherIsBetter = FALSE)

  saveSection(2, result2, list(
    summary = list(table = headline, caption = paste(
      "Sojourn recovery: Viterbi ARI and TV distance between fitted and true dwell pmfs (mean over states),",
      "mean (sd) over replicates. THMM rows show the geometric dwell a THMM implies for the same data.")),
    dwell = list(table = dwellTable, caption = "Sojourn recovery: dwell-time parameter estimates of matched states."),
    emission = list(table = emissionTable, caption = "Sojourn recovery: emission parameter estimates of matched states."),
    pairedTV = list(table = pairedTv, caption = paste(
      "Sojourn recovery: paired reduction in dwell TV distance (THMM minus THSMM) on the same data sets."))
  ))
  printTable("Sojourn recovery: mean (sd) over replicates", headline)
  printTable("Dwell parameters, true family (bias / RMSE)",
             dwellTable[dwellTable$method == "THSMM, true family" & dwellTable$parameter != "TV dwell" &
                          !grepl("curves", dwellTable$scenario), ])

  pdf(file.path(outputDir, "section2.pdf"), width = 11, height = 4)
  for (name in allProcesses) {
    for (scenarioName in names(section2Scenarios)) {
      s <- section2Scenarios[[scenarioName]]
      if (!name %in% s$processes || !is.na(s$series)) next
      truth <- makeTruth(processes[[name]], s)
      p <- result2$pmf[result2$pmf$process == name & result2$pmf$scenario == scenarioName, ]
      par(mfrow = c(1, truth$nStates), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0), las = 1)
      for (i in seq_len(truth$nStates)) {
        shown <- setdiff(names(s$params[[i]]), "dwellSize")
        label <- if (s$type[i] == "np") sprintf("State %d: bimodal at %s", i, paste(which(s$params[[i]] > 0), collapse = " and "))
                 else sprintf("State %d: %s", i, paste(shown, signif(s$params[[i]][shown], 3), sep = " = ", collapse = ", "))
        plotDwellPanel(p, truth, i, label, "THSMM, true family")
        if (i == 1) methodLegend("topright", intersect(methodOrder, unique(p$method)), cex = 0.7)
      }
      mtext(sprintf("Section 2, %s dwell times: %s", tolower(scenarioName), processes[[name]]$label), outer = TRUE, font = 2)
    }
  }
  # Sample-size series: RMSE of the mean dwell and mean TV distance vs number of curves.
  series <- names(section2Scenarios)[vapply(section2Scenarios, function(s) identical(s$series, "Sample size"), logical(1))]
  seriesRows <- rows[rows$process == "BM" & rows$scenario %in% c("Overdispersed", series), ]
  seriesRows$n <- vapply(seriesRows$scenario, function(s) section2Scenarios[[s]]$nCurves, numeric(1))
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0), las = 1)
  for (metric in c("TV dwell", "mean dwell")) {
    g <- seriesRows[seriesRows$parameter == metric, ]
    g$value <- if (metric == "TV dwell") g$estimate else (g$estimate - g$truth)^2
    a <- aggregate(value ~ n + method, g, mean)
    if (metric == "mean dwell") a$value <- sqrt(a$value)
    methods <- intersect(methodOrder, unique(a$method))
    plot(range(a$n), range(a$value), type = "n", log = "x", xlab = "number of curves",
         ylab = if (metric == "TV dwell") "mean TV distance" else "RMSE of mean dwell",
         main = if (metric == "TV dwell") "Dwell pmf" else "Mean dwell", font.main = 1)
    for (m in methods) with(a[a$method == m, ][order(a$n[a$method == m]), ], lines(n, value, type = "b", pch = 19, col = methodColours[[m]]))
    if (metric == "TV dwell") methodLegend("topright", methods, truth = FALSE, cex = 0.75)
  }
  mtext("Section 2, sample-size series: overdispersed dwell times, Brownian motion", outer = TRUE, font = 2)
  invisible(dev.off())
}


# ===========================================================================
# Section 3: Viterbi ARI, THMM vs THSMM, non-geometric dwell times
# ===========================================================================

if (3 %in% sections) {
  cat("\nSection 3: Viterbi ARI under non-geometric dwell times\n")
  methods3 <- c("THMM (funHMM)", "THSMM, geometric", "THSMM, npg(12)", "THSMM, true family")
  result3 <- runTasks(makeTasks(section3Scenarios), replicateFunction(3, section3Scenarios, methods3), 3)
  rows <- result3$rows

  ari <- rows[rows$parameter == "ARI", ]
  ariTable <- do.call(rbind, lapply(split(ari, ari[, c("process", "scenario")], drop = TRUE), function(g) {
    methods <- intersect(methodOrder, unique(g$method))
    cells <- vapply(methods, function(m) meanSd(g$estimate[g$method == m]), character(1))
    failed <- rows[rows$process == g$process[1] & rows$scenario == g$scenario[1] & rows$parameter == "failed", ]
    cbind(g[1, c("process", "scenario")], as.data.frame(t(cells), check.names = FALSE),
          failures = if (nrow(failed)) paste(table(failed$method), names(table(failed$method)), collapse = "; ") else "")
  }))
  rownames(ariTable) <- NULL
  ariTable <- ariTable[order(match(ariTable$process, allProcesses), ariTable$scenario), ]
  paired <- pairedTable(rows, "ARI")

  saveSection(3, result3, list(
    ari = list(table = ariTable, caption = paste(
      "Viterbi ARI under non-geometric dwell times, mean (sd) over replicates. THSMM, geometric is the",
      "THMM fitted by the THSMM code; failed fits are excluded.")),
    paired = list(table = paired, caption = paste(
      "Paired ARI differences, THSMM minus THMM, on the same data sets: mean, standard error, share of",
      "replicates where the THSMM decodes better, and Wilcoxon signed-rank p-value."))
  ))
  printTable("Viterbi ARI: mean (sd) over replicates", ariTable)
  printTable("Paired ARI gains of the THSMM", paired)

  pdf(file.path(outputDir, "section3.pdf"), width = 11, height = 4)
  metricBoxplots(rows, "ARI", section3Scenarios, "Viterbi ARI", "Section 3, Viterbi ARI", ylim = c(0, 1))

  # Sensitivity series (bimodal truth, BM): mean ARI +- 2 se.
  bimodalRows <- ari[ari$process == "BM" & ari$scenario %in% names(section3Scenarios)[
    vapply(section3Scenarios, function(s) !is.na(s$series) || identical(s$type, bimodalDwell$type), logical(1))], ]
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0), las = 1)
  for (seriesName in c("Separation", "Sequence length")) {
    members <- names(section3Scenarios)[vapply(section3Scenarios, function(s) {
      identical(s$type, bimodalDwell$type) && (identical(s$series, seriesName) ||
        (s$level == "low" && s$nCurves == round(500 * curveScale)))
    }, logical(1))]
    value <- vapply(members, function(s) {
      s <- section3Scenarios[[s]]
      if (seriesName == "Separation") abs(processes$BM$params[[s$level]][1, 1]) else s$nCurves
    }, numeric(1))
    g <- bimodalRows[bimodalRows$scenario %in% members, ]
    g$value <- value[g$scenario]
    a <- do.call(rbind, lapply(split(g, g[, c("value", "method")], drop = TRUE), function(d) {
      data.frame(value = d$value[1], method = d$method[1], mean = mean(d$estimate), se = sd(d$estimate) / sqrt(nrow(d)))
    }))
    methods <- intersect(methodOrder, unique(a$method))
    plot(range(a$value), c(0, 1), type = "n", log = if (seriesName == "Sequence length") "x" else "",
         xlab = if (seriesName == "Separation") "drift magnitude |c| (states -c, 0, c)" else "number of curves",
         ylab = "Viterbi ARI", main = seriesName, font.main = 1)
    for (m in methods) {
      d <- a[a$method == m, ][order(a$value[a$method == m]), ]
      arrows(d$value, d$mean - 2 * d$se, d$value, d$mean + 2 * d$se, angle = 90, code = 3, length = 0.03, col = methodColours[[m]])
      lines(d$value, d$mean, type = "b", pch = 19, col = methodColours[[m]])
    }
    if (seriesName == "Separation") methodLegend("bottomright", methods, truth = FALSE, cex = 0.75)
  }
  mtext("Section 3, bimodal dwell times, Brownian motion: sensitivity", outer = TRUE, font = 2)
  invisible(dev.off())
}

cat(sprintf("\nResults written to %s\n", outputDir))
