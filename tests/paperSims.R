# Simulation studies for the THSMM paper.
#
# Run from the package root:
#   Rscript tests/paperSims.R            all sections (0-3)
#   Rscript tests/paperSims.R 1 3        selected sections (0: example curves only, seconds)
#   Rscript tests/paperSims.R --quick    2 replicates, fewer curves: a smoke test
# Environment variables: PAPER_SIMS_REPS (replicates per setting, default 50)
# and PAPER_SIMS_CORES (parallel workers, default detectCores() - 1).
#
# Every section is run for the three emission settings of Kashlak, Loliencar
# and Heo (2023, Sec. 5), with K = 3 states:
#   BM    Brownian motion with linear drift c           ("bmwd", hurst 1/2)
#   OU    Ornstein-Uhlenbeck dO = (b0 - b1 O) dtau + dW  ("ou", b0 = c mu, b1 = c)
#   fBM   fractional BM, hurst 0.75, linear drift c      ("bmwd", hurst known)
# Curves are observed at 101 equispaced points of [0, 1] with O(0) = 0. The
# THMM package funHMM puts its len = 100 points at (1:len) / len, so dropping
# the first column gives funHMM the data layout its own rbm() / rou() produce
# for the same process.
#
# Every THMM is fitted with funHMM::thmm() (Baum-Welch in C) with sigma =
# thmmSigma = 1, the unit diffusion the data are simulated with and the
# package's bmwd / ou functionals assume, so THMM and THSMM emissions are on the
# same scale. It starts from the package's k-means start for a geometric dwell
# (emission parameters, transition matrix, initial distribution), so THMM and
# THSMM share the emission start. THSMMs are fitted with this package.
#
# Section 0 -- figure of example emission curves for each process at low and
# medium separation (example-curves.pdf).
#
# Section 1 -- THMM recovery. Data from a THMM (geometric dwell times), fitted by
#   * THMM (funHMM)        funHMM::thmm()
#   * THSMM, npg(10)       nonparametric dwell: should find the geometric shape
#   from a common emission start. Scored on emission and dwell parameters,
#   transition matrix, dwell pmfs, Viterbi ARI against the truth, and agreement
#   of the THSMM and THMM Viterbi paths.
#
# Section 2 -- sojourn recovery. Data from THSMMs whose states share one dwell
# family (geometric, Poisson, shifted Poisson, negative binomial, Weibull),
# fitted with the true family (shifts profiled by fitHsmmShiftGrid) and with
# npg(10). Bias / sd / RMSE of every dwell and emission parameter, the mean
# dwell, and the total-variation distance between fitted and true dwell pmfs.
#
# Section 3 -- Viterbi ARI, THMM vs THSMM, when dwell times are not geometric
# (shifted Poisson, peaked negative binomial, Weibull with increasing hazard,
# and a truth mixing the three families), with emissions at low separation so
# that the dwell structure matters for decoding. Methods: THMM (funHMM),
# THSMM npg(10), THSMM with the true family.
#
# Output in tests/output/paperSims/: per-replicate results (sectionN.rds and
# sectionN-replicates.csv), summary tables (.csv and LaTeX .tex) and figures
# (sectionN.pdf).

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

if (!requireNamespace("funHMM", quietly = TRUE)) stop("funHMM is not installed")

library(parallel)

arguments <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% arguments
sections <- as.integer(setdiff(arguments, "--quick"))
if (length(sections) == 0L) sections <- 0:3

nReps <- if (quick) 2L else as.integer(Sys.getenv("PAPER_SIMS_REPS", "50"))
nCores <- as.integer(Sys.getenv("PAPER_SIMS_CORES", max(1L, detectCores() - 1L)))

outputDir <- file.path("tests", "output", "paperSims")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

K <- 3
grid <- seq(0, 1, length.out = 101)
nCurves <- if (quick) c(200, 300, 200) else c(500, 1000, 500)   # per section
npgSize <- 10                     # npg: 9 free head probabilities + geometric tail
parametricDwellSize <- 30         # exact dwell support before the geometric tail
shiftCandidates <- if (quick) 2:4 else 1:4
maxDwell <- 300                   # support for mean dwells and TV distances
maxPlot <- 20                     # dwell support shown in the figures
emTol <- 1e-8
emMaxIter <- 2000
thmmTol <- 1e-8                   # funHMM::thmm() Baum-Welch
thmmMaxIter <- 1000
thmmSigma <- 1                    # unit diffusion, as simulated and as the package's functionals

ouParams <- function(mu, rate) cbind(b0 = rate * mu, b1 = rate)

# Emission parameters at two separations. Single-curve classification accuracy
# at the true parameters is roughly 0.9+ at "medium" and 0.78 at "low" for each
# process, so with "low" the dwell structure carries real information. States
# are ordered by drift / mean level, as initializeHsmmEM() orders its starts.
processes <- list(
  BM = list(
    label = "Brownian motion with drift", emissionType = "bmwd", hurst = 0.5,
    control = list(grid = grid),
    params = list(medium = cbind(drift = c(-3, 0, 3)), low = cbind(drift = c(-2, 0, 2)))
  ),
  OU = list(
    label = "Ornstein-Uhlenbeck", emissionType = "ou", hurst = 0.5,
    control = list(grid = grid),
    params = list(medium = ouParams(c(-1, 0, 1), 3), low = ouParams(c(-1, 0, 1), 2))
  ),
  fBM = list(
    label = "Fractional BM with drift (hurst 0.75)", emissionType = "bmwd", hurst = 0.75,
    control = list(grid = grid, hurst = 0.75),
    params = list(medium = cbind(drift = c(-3, 0, 3)), low = cbind(drift = c(-2, 0, 2)))
  )
)

uniformOmega <- matrix(1 / (K - 1), K, K)
diag(uniformOmega) <- 0

# Section 1: THMM with leave probabilities 0.4, 0.25, 0.15 (mean dwells 2.5, 4,
# 6.7) and uniform jumps, i.e. A = diag(1 - p) + diag(p) omega.
thmmLeave <- c(0.4, 0.25, 0.15)

withSize <- function(...) c(..., dwellSize = parametricDwellSize)

# Section 2: one family per scenario, distinct parameters per state.
sojournScenarios <- list(
  Geometric = list(
    type = rep("geo", K),
    params = list(c(prob = 0.5), c(prob = 0.3), c(prob = 0.15))                       # means 2, 3.3, 6.7
  ),
  Poisson = list(
    type = rep("pois", K),
    params = list(withSize(lambda = 2), withSize(lambda = 5), withSize(lambda = 9))    # 1 + Pois: means 3, 6, 10
  ),
  `Shifted Poisson` = list(
    type = rep("shiftp", K),
    params = list(withSize(lambda = 1, shift = 3), withSize(lambda = 3, shift = 2),
                  withSize(lambda = 2, shift = 4))                                     # means 4, 5, 6
  ),
  `Negative binomial` = list(
    type = rep("nb", K),
    params = list(withSize(size = 2, prob = 0.4), withSize(size = 5, prob = 0.5),
                  withSize(size = 10, prob = 0.55))                                    # 1 + NB: means 4, 6, 9.2
  ),
  Weibull = list(
    type = rep("weibull", K),
    params = list(withSize(shape = 1.5, scale = 3), withSize(shape = 2.5, scale = 6),
                  withSize(shape = 4, scale = 9))                                      # means 3.2, 5.8, 8.7
  )
)

# Section 3: non-geometric dwell times with a typical length (mode > 1).
nonGeometricScenarios <- list(
  `Shifted Poisson` = list(
    type = rep("shiftp", K),
    params = list(withSize(lambda = 1, shift = 3), withSize(lambda = 2, shift = 3),
                  withSize(lambda = 1, shift = 4))                                     # means 4, 5, 5
  ),
  `Peaked negative binomial` = list(
    type = rep("nb", K),
    params = list(withSize(size = 20, prob = 20 / 23), withSize(size = 20, prob = 20 / 24),
                  withSize(size = 20, prob = 20 / 25))                                 # means 4, 5, 6
  ),
  `Weibull (increasing hazard)` = list(
    type = rep("weibull", K),
    params = list(withSize(shape = 3, scale = 4), withSize(shape = 3, scale = 5),
                  withSize(shape = 3, scale = 6))                                      # means 4.1, 5, 5.9
  ),
  `Mixed families` = list(
    type = c("shiftp", "nb", "weibull"),
    params = list(withSize(lambda = 1, shift = 3), withSize(size = 20, prob = 20 / 24),
                  withSize(shape = 3, scale = 6))
  )
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

adjustedRand <- function(a, b) mclust::adjustedRandIndex(a, b)

permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  smaller <- permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))))
}
allPermutations <- permutations(K)

# match[i] = estimated state that best matches true state i (maximizes matched counts).
matchStates <- function(estimated, truth) {
  tab <- table(factor(estimated, seq_len(K)), factor(truth, seq_len(K)))
  scores <- apply(allPermutations, 1, function(p) sum(tab[cbind(p, seq_len(K))]))
  allPermutations[which.max(scores), ]
}

seedFor <- function(section, processIndex, scenarioIndex, replicate) {
  as.integer(1e6 * section + 1e4 * processIndex + 1e2 * scenarioIndex + replicate)
}

makeTruth <- function(process, sojournType, sojournParameters, level) {
  buildHsmm(uniformOmega, sojournParameters, sojournType, process$params[[level]], process$emissionType,
            delta = rep(1 / K, K), stationary = FALSE, emissionControl = process$control)
}

thmmTruth <- function(process) {
  makeTruth(process, rep("geo", K), lapply(thmmLeave, function(p) c(prob = p)), "medium")
}

# Starting model from the data alone; deterministic given the seed, so every
# dwell family shares the same emission start.
startModel <- function(x, process, sojournType, seed) {
  sojournType <- rep_len(sojournType, K)
  dwellSize <- ifelse(sojournType == "npg", npgSize, parametricDwellSize)
  initializeHsmmEM(x, K, sojournType, process$emissionType, process$control,
                   dwellSize = dwellSize, stationary = FALSE, seed = seed)$model
}

fitThsmm <- function(start, x) {
  if (any(start$sojournType %in% c("shiftp", "shiftnb"))) {
    # Shifts are integers held fixed by EM: profile them over a grid.
    shiftValues <- lapply(start$sojournType, function(type) if (type %in% c("shiftp", "shiftnb")) shiftCandidates)
    return(suppressWarnings(fitHsmmShiftGrid(start, x, shiftValues, stationary = FALSE, fitFunction = fitHsmmEM,
                                             tol = emTol, maxIter = emMaxIter)))
  }
  suppressWarnings(fitHsmmEM(start, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
}

# THMM transition matrix of a geometric-dwell model: A = diag(1 - p) + diag(p) omega.
thmmMatrix <- function(model) {
  leave <- vapply(model$sojournParameters, function(s) s[["prob"]], numeric(1))
  diag(1 - leave) + leave * model$tpm
}

# Package emission parameters <-> funHMM's `par`: drifts for bmwd; for ou the
# package's (b0, b1) = (rate * mean, rate) against funHMM's (mean, rate). A zero
# rate (the OU boundary) is floored so that mean * rate still returns b0.
toFunhmmPar <- function(params, emissionType) {
  if (emissionType != "ou") return(params[, 1])
  rate <- pmax(params[, 2], 1e-8)
  cbind(mean = params[, 1] / rate, rate = rate)
}
fromFunhmmPar <- function(par, emissionType) {
  if (emissionType != "ou") return(matrix(par, ncol = 1, dimnames = list(NULL, "drift")))
  # funHMM reports mean = NA when a fitted rate is 0, so b0 is then unknown.
  cbind(intercept = par[, "mean"] * par[, "rate"], rate = par[, "rate"])
}

# funHMM::thmm() Baum-Welch (THMM) from a package geometric-dwell starting model.
fitFunhmm <- function(x, start, process) {
  fit <- tryCatch(suppressWarnings(funHMM::thmm(
    x[, -1], K, type = process$emissionType, hurst = process$hurst, sigma = thmmSigma,
    par = toFunhmmPar(start$emissionParameters, process$emissionType),
    trans = thmmMatrix(start), init.prob = start$delta,
    max.iter = thmmMaxIter, min.iter = 10, tol = thmmTol
  )), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(states = fit$states, A = fit$trans, emissionParameters = fromFunhmmPar(fit$par, process$emissionType))
}

# Comparable emission summaries: drift, or (b0, b1, mean level mu) for OU.
emissionSummary <- function(params, emissionType) {
  if (emissionType == "ou") return(cbind(b0 = params[, 1], b1 = params[, 2], mu = params[, 1] / params[, 2]))
  cbind(drift = params[, 1])
}

# A fitted method in a common form: decoded states, emission parameters, and
# per-(estimated)-state dwell pmf / dwell family / dwell parameters.
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
  list(
    states = fit$states,
    emissionParameters = fit$emissionParameters,
    dwellPmf = function(j) dgeom(seq_len(maxDwell) - 1, 1 - fit$A[j, j]),
    dwellType = function(j) "geo",
    dwellParams = function(j) c(prob = 1 - fit$A[j, j]),
    A = fit$A
  )
}

resultRow <- function(state, parameter, truth, estimate) {
  data.frame(state = state, parameter = parameter, truth = truth, estimate = estimate)
}

# Score one fitted method against the truth: ARI, emission and dwell parameters
# of matched states, mean dwell, TV distance of the dwell pmf, transition
# matrix error. Also returns the matched dwell pmfs on 1:maxPlot for figures.
scoreMethod <- function(method, fitted, truth, trueStates) {
  match <- matchStates(fitted$states, trueStates)
  rows <- list(resultRow(NA, "ARI", 1, adjustedRand(fitted$states, trueStates)))

  trueEmission <- emissionSummary(truth$emissionParameters, truth$emissionType)
  estEmission <- emissionSummary(fitted$emissionParameters, truth$emissionType)[match, , drop = FALSE]
  dwellTimes <- seq_len(maxDwell)
  pmfRows <- list()

  for (i in seq_len(K)) {
    j <- match[i]
    for (p in colnames(trueEmission)) rows[[length(rows) + 1L]] <- resultRow(i, p, trueEmission[i, p], estEmission[i, p])

    truePmf <- getModelDwellPmf(truth, i, maxDwell)
    estPmf <- fitted$dwellPmf(j)
    rows[[length(rows) + 1L]] <- resultRow(i, "mean dwell", sum(dwellTimes * truePmf), sum(dwellTimes * estPmf))
    rows[[length(rows) + 1L]] <- resultRow(i, "TV dwell", 0, 0.5 * sum(abs(truePmf - estPmf)))

    # Dwell parameters, when the fitted family is the true one (dwellSize is fixed).
    if (fitted$dwellType(j) == truth$sojournType[i]) {
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

# Run fitted methods (a named list of functions returning the common form, or
# NULL on failure) and score each.
scoreMethods <- function(methods, truth, sim) {
  fitted <- list()
  scores <- list()
  seconds <- numeric()
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
  list(fitted = fitted,
       rows = do.call(rbind, lapply(scores, `[[`, "rows")),
       pmf = do.call(rbind, lapply(scores, `[[`, "pmf")))
}

# Run replicate tasks in parallel. Each task returns list(rows, pmf); the
# identifying columns of the task are prepended.
runTasks <- function(tasks, replicateFunction, section) {
  started <- Sys.time()
  outputs <- mclapply(seq_len(nrow(tasks)), function(t) {
    task <- tasks[t, , drop = FALSE]
    tryCatch(replicateFunction(task), error = function(e) list(error = conditionMessage(e)))
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

# Bias / sd / RMSE of every (method, state, parameter) across replicates.
summarizeEstimates <- function(rows, by) {
  rows <- rows[!is.na(rows$truth) & rows$parameter != "seconds", ]
  groups <- split(rows, rows[, by], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    error <- g$estimate - g$truth
    cbind(g[1, by, drop = FALSE],
          data.frame(truth = g$truth[1], mean = mean(g$estimate, na.rm = TRUE), bias = mean(error, na.rm = TRUE),
                     sd = sd(g$estimate, na.rm = TRUE), rmse = sqrt(mean(error^2, na.rm = TRUE)),
                     n = sum(!is.na(g$estimate))))
  }))
  rownames(out) <- NULL
  out
}

meanSd <- function(v, digits = 3) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_character_)
  sprintf("%.*f (%.*f)", digits, mean(v), digits, sd(v))
}

writeLatex <- function(table, file, caption, digits = 3) {
  escape <- function(s) gsub("([&%#_$])", "\\\\\\1", s)
  cells <- vapply(table, function(column) {
    if (is.numeric(column)) formatC(column, digits = digits, format = "f") else escape(as.character(column))
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
    writeLatex(tables[[name]]$table, file.path(outputDir, sprintf("section%d-%s.tex", section, name)),
               tables[[name]]$caption)
  }
}

methodColours <- c("THMM (funHMM)" = "#4d4d4d", "THSMM, npg(10)" = "#b2182b", "THSMM, true family" = "#1a9850")

# Dwell pmfs: truth (grey bars), fitted median per method (lines) and the
# 10-90% band of the first method listed in `bands`.
plotDwellPanel <- function(pmf, truth, i, main, bands) {
  truePmf <- getModelDwellPmf(truth, i, maxPlot)
  methods <- intersect(names(methodColours), unique(pmf$method))
  quantiles <- lapply(methods, function(m) {
    d <- pmf[pmf$method == m & pmf$state == i, ]
    sapply(split(d$prob, d$dwell), quantile, c(0.1, 0.5, 0.9), na.rm = TRUE)
  })
  names(quantiles) <- methods
  upper <- max(truePmf, unlist(lapply(quantiles, function(q) q[3, ])))
  plot(seq_len(maxPlot), truePmf, type = "h", lwd = 6, col = "grey75", ylim = c(0, upper * 1.05),
       xlab = "dwell time (curves)", ylab = "probability", main = main)
  for (m in intersect(bands, methods)) {
    polygon(c(seq_len(maxPlot), rev(seq_len(maxPlot))), c(quantiles[[m]][1, ], rev(quantiles[[m]][3, ])),
            col = adjustcolor(methodColours[[m]], 0.2), border = NA)
  }
  for (m in methods) lines(seq_len(maxPlot), quantiles[[m]][2, ], type = "b", pch = 19, cex = 0.5, col = methodColours[[m]])
}

dwellLegend <- function(methods) {
  legend("topright", c("truth", methods), col = c("grey75", methodColours[methods]),
         lwd = c(6, rep(1.5, length(methods))), bty = "n", cex = 0.75)
}

cat(sprintf("paperSims: sections %s, %d replicates, %d cores%s\n", paste(sections, collapse = ", "),
            nReps, nCores, if (quick) " (quick mode)" else ""))


# ===========================================================================
# Section 0: example emission curves at low and medium separation
# ===========================================================================

# Rows: low / medium separation; columns: BM, OU, fBM. A few simulated curves
# per state (thin lines) and each state's mean function E O(tau) (thick):
# c tau for BM and fBM, mu (1 - exp(-c tau)) for OU started at 0.
if (0 %in% sections) {
  cat("\nSection 0: example emission curves\n")

  curvesPerState <- 8
  stateColours <- c("#0072B2", "#D55E00", "#009E73")   # Okabe-Ito, colour-blind safe
  levels <- c(low = "Low separation", medium = "Medium separation")
  states <- rep(seq_len(K), each = curvesPerState)

  meanFunction <- function(process, params, i) {
    if (process$emissionType == "ou") return(ouMeanLevel(params)[i, 1] * (1 - exp(-params[i, 2] * grid)))
    params[i, 1] * grid
  }

  parameterText <- function(process, params) {
    if (process$emissionType == "ou") {
      return(sprintf("mu = (%s), c = %g", paste(ouMeanLevel(params)[, 1], collapse = ", "), params[1, 2]))
    }
    sprintf("c = (%s)", paste(params[, 1], collapse = ", "))
  }

  set.seed(2026)
  examples <- lapply(processes, function(process) lapply(names(levels), function(level) {
    family <- getEmissionFamily(process$emissionType)
    family$simulate(states, process$params[[level]], process$control)
  }))

  pdf(file.path(outputDir, "example-curves.pdf"), width = 7.5, height = 4.8, pointsize = 10)
  par(mfrow = c(2, length(processes)), mar = c(3.2, 3.2, 1.6, 0.6), oma = c(2.4, 2, 0.4, 0),
      mgp = c(2, 0.6, 0), tcl = -0.3, las = 1)

  for (row in seq_along(levels)) {
    for (name in names(processes)) {
      process <- processes[[name]]
      params <- process$params[[names(levels)[row]]]
      x <- examples[[name]][[row]]
      yRange <- range(unlist(lapply(examples[[name]], range)))   # shared within a column

      matplot(grid, t(x), type = "l", lty = 1, lwd = 0.6, col = adjustcolor(stateColours[states], 0.55),
              ylim = yRange, xlab = expression(tau), ylab = if (name == names(processes)[1]) expression(O(tau)) else "",
              main = if (row == 1) c(BM = "Brownian motion with drift", OU = "Ornstein-Uhlenbeck",
                                     fBM = "Fractional BM with drift (H = 0.75)")[[name]] else "",
              font.main = 1, cex.main = 1)
      for (i in seq_len(K)) lines(grid, meanFunction(process, params, i), col = stateColours[i], lwd = 2.2)
      legend("topleft", parameterText(process, params), bty = "n", cex = 0.8, inset = c(-0.04, -0.02))
    }
    mtext(levels[[row]], side = 2, outer = TRUE, line = 0.6, at = c(0.76, 0.27)[row], las = 0)
  }

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", c(paste("State", seq_len(K)), "Mean function"), col = c(stateColours, "grey30"),
         lwd = c(rep(1.2, K), 2.2), horiz = TRUE, bty = "n", cex = 0.9)
  invisible(dev.off())
  cat(sprintf("  Wrote %s\n", file.path(outputDir, "example-curves.pdf")))
}


# ===========================================================================
# Section 1: can the THSMM recover the THMM?
# ===========================================================================

if (1 %in% sections) {
  cat("\nSection 1: THMM recovery\n")

  section1Replicate <- function(task) {
    process <- processes[[task$process]]
    seed <- seedFor(1, match(task$process, names(processes)), 1, task$replicate)
    truth <- thmmTruth(process)
    sim <- simulateHsmm(truth, nCurves[1], initialType = "aggregateStart", seed = seed)
    geoStart <- startModel(sim$x, process, "geo", seed)
    npgStart <- startModel(sim$x, process, "npg", seed)

    scored <- scoreMethods(list(
      "THMM (funHMM)" = function() { fit <- fitFunhmm(sim$x, geoStart, process); if (!is.null(fit)) fromThmmFit(fit) },
      "THSMM, npg(10)" = function() fromPackageFit(fitThsmm(npgStart, sim$x), sim$x)
    ), truth, sim)

    # Agreement of the THSMM decoding with the THMM's.
    reference <- scored$fitted[["THMM (funHMM)"]]
    agreement <- if (!is.null(reference)) do.call(rbind, lapply(setdiff(names(scored$fitted), "THMM (funHMM)"), function(m) {
      data.frame(method = m, state = NA, parameter = "ARI vs THMM", truth = NA,
                 estimate = adjustedRand(scored$fitted[[m]]$states, reference$states))
    }))
    list(rows = rbind(scored$rows, agreement), pmf = scored$pmf)
  }

  tasks <- expand.grid(replicate = seq_len(nReps), process = names(processes), stringsAsFactors = FALSE)[, 2:1]
  result1 <- runTasks(tasks, section1Replicate, 1)
  rows <- result1$rows

  estimates <- summarizeEstimates(rows, c("process", "method", "state", "parameter"))

  headline <- do.call(rbind, lapply(split(rows, rows[, c("process", "method")], drop = TRUE, lex.order = TRUE), function(g) {
    pick <- function(p) g$estimate[g$parameter == p]
    data.frame(process = g$process[1], method = g$method[1],
               ARI = meanSd(pick("ARI")),
               `ARI vs THMM` = meanSd(pick("ARI vs THMM")),
               `mean TV dwell` = meanSd(pick("TV dwell")),
               `max |A error|` = meanSd(pick("max |A error|")),
               failures = sum(g$parameter == "failed"),
               seconds = round(mean(pick("seconds")), 1), check.names = FALSE)
  }))
  rownames(headline) <- NULL

  saveSection(1, result1, list(
    summary = list(table = headline, caption = paste(
      "THMM recovery: Viterbi ARI against the truth and against the THMM (funHMM), total-variation distance",
      "between fitted and true (geometric) dwell pmfs, and max abs error of the THMM transition matrix;",
      "mean (sd) over replicates.")),
    estimates = list(table = estimates, caption = "THMM recovery: parameter estimates of matched states.")
  ))

  cat("\n")
  print(headline, row.names = FALSE)
  cat("\n  Emission and dwell parameters (bias / RMSE):\n")
  print(format(estimates[!estimates$parameter %in% c("ARI", "TV dwell", "max |A error|"), ], digits = 3), row.names = FALSE)

  pdf(file.path(outputDir, "section1.pdf"), width = 11, height = 7)
  for (name in names(processes)) {
    process <- processes[[name]]
    truth <- thmmTruth(process)
    r <- rows[rows$process == name, ]
    p <- result1$pmf[result1$pmf$process == name, ]
    methods <- intersect(names(methodColours), unique(r$method))
    par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

    ari <- r[r$parameter == "ARI", ]
    boxplot(estimate ~ factor(method, methods), data = ari, col = methodColours[methods], names = seq_along(methods),
            xlab = "method", ylab = "ARI", main = "Viterbi ARI vs truth")
    legend("bottomleft", paste(seq_along(methods), methods), bty = "n", cex = 0.7)

    firstParameter <- colnames(emissionSummary(truth$emissionParameters, truth$emissionType))[if (truth$emissionType == "ou") 3 else 1]
    emission <- r[r$parameter == firstParameter, ]
    boxplot(estimate - truth ~ factor(method, methods) + state, data = emission, col = methodColours[methods],
            xlab = "method . state", ylab = paste("error in", firstParameter), main = paste("Emission:", firstParameter), las = 2,
            cex.axis = 0.6)
    abline(h = 0, lty = 2)

    stay <- r[r$parameter == "prob", ]
    boxplot(estimate - truth ~ factor(method, methods) + state, data = stay, col = methodColours[methods],
            xlab = "method . state", ylab = "error in leave probability", main = "Dwell: leave probability 1 - A_ii",
            las = 2, cex.axis = 0.6)
    abline(h = 0, lty = 2)

    for (i in seq_len(K)) {
      plotDwellPanel(p, truth, i, sprintf("Dwell pmf, state %d (geometric, p = %.2f)", i, thmmLeave[i]), "THSMM, npg(10)")
      if (i == 1) dwellLegend(methods)
    }
    mtext(paste("Section 1, THMM recovery:", process$label), outer = TRUE, font = 2)
  }
  invisible(dev.off())
}


# ===========================================================================
# Section 2: recovery of sojourn distributions
# ===========================================================================

if (2 %in% sections) {
  cat("\nSection 2: sojourn recovery\n")

  section2Replicate <- function(task) {
    process <- processes[[task$process]]
    scenario <- sojournScenarios[[task$scenario]]
    seed <- seedFor(2, match(task$process, names(processes)), match(task$scenario, names(sojournScenarios)), task$replicate)
    truth <- makeTruth(process, scenario$type, scenario$params, "medium")
    sim <- simulateHsmm(truth, nCurves[2], initialType = "aggregateStart", seed = seed)

    scored <- scoreMethods(list(
      "THSMM, true family" = function() fromPackageFit(fitThsmm(startModel(sim$x, process, scenario$type, seed), sim$x), sim$x),
      "THSMM, npg(10)" = function() fromPackageFit(fitThsmm(startModel(sim$x, process, "npg", seed), sim$x), sim$x)
    ), truth, sim)
    scored[c("rows", "pmf")]
  }

  tasks <- expand.grid(replicate = seq_len(nReps), scenario = names(sojournScenarios), process = names(processes),
                       stringsAsFactors = FALSE)[, 3:1]
  result2 <- runTasks(tasks, section2Replicate, 2)
  rows <- result2$rows

  estimates <- summarizeEstimates(rows, c("process", "scenario", "method", "state", "parameter"))
  estimates <- estimates[estimates$parameter != "ARI", ]
  dwellParameters <- c("prob", "lambda", "shift", "size", "shape", "scale", "mean dwell", "TV dwell")
  dwellTable <- estimates[estimates$parameter %in% dwellParameters, ]
  emissionTable <- estimates[!estimates$parameter %in% dwellParameters, ]

  headline <- do.call(rbind, lapply(split(rows, rows[, c("process", "scenario", "method")], drop = TRUE, lex.order = TRUE), function(g) {
    pick <- function(p) g$estimate[g$parameter == p]
    meanDwell <- g[g$parameter == "mean dwell", ]
    data.frame(process = g$process[1], scenario = g$scenario[1], method = g$method[1],
               ARI = meanSd(pick("ARI")),
               `mean TV dwell` = meanSd(pick("TV dwell")),
               `mean dwell RMSE` = round(sqrt(mean((meanDwell$estimate - meanDwell$truth)^2)), 3),
               failures = sum(g$parameter == "failed"), check.names = FALSE)
  }))
  rownames(headline) <- NULL

  saveSection(2, result2, list(
    summary = list(table = headline, caption = paste(
      "Sojourn recovery: Viterbi ARI, total-variation distance between fitted and true dwell pmfs",
      "(mean (sd) over replicates and states), and RMSE of the state mean dwells.")),
    dwell = list(table = dwellTable, caption = "Sojourn recovery: dwell-time parameter estimates of matched states."),
    emission = list(table = emissionTable, caption = "Sojourn recovery: emission parameter estimates of matched states.")
  ))

  cat("\n")
  print(headline, row.names = FALSE)
  cat("\n  Dwell parameters, true family (bias / RMSE):\n")
  print(format(dwellTable[dwellTable$method == "THSMM, true family" & !dwellTable$parameter %in% c("TV dwell"), ],
               digits = 3), row.names = FALSE)

  pdf(file.path(outputDir, "section2.pdf"), width = 11, height = 4)
  for (name in names(processes)) {
    process <- processes[[name]]
    for (scenarioName in names(sojournScenarios)) {
      scenario <- sojournScenarios[[scenarioName]]
      truth <- makeTruth(process, scenario$type, scenario$params, "medium")
      p <- result2$pmf[result2$pmf$process == name & result2$pmf$scenario == scenarioName, ]
      par(mfrow = c(1, K), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
      for (i in seq_len(K)) {
        shown <- setdiff(names(scenario$params[[i]]), "dwellSize")
        plotDwellPanel(p, truth, i, sprintf("State %d: %s", i, paste(shown, signif(scenario$params[[i]][shown], 3),
                                                                     sep = " = ", collapse = ", ")),
                       "THSMM, true family")
        if (i == 1) dwellLegend(c("THSMM, true family", "THSMM, npg(10)"))
      }
      mtext(sprintf("Section 2, %s dwell times: %s", scenarioName, process$label), outer = TRUE, font = 2)
    }
  }
  invisible(dev.off())
}


# ===========================================================================
# Section 3: Viterbi ARI, THMM vs THSMM, non-geometric dwell times
# ===========================================================================

if (3 %in% sections) {
  cat("\nSection 3: Viterbi ARI under non-geometric dwell times\n")

  section3Replicate <- function(task) {
    process <- processes[[task$process]]
    scenario <- nonGeometricScenarios[[task$scenario]]
    seed <- seedFor(3, match(task$process, names(processes)), match(task$scenario, names(nonGeometricScenarios)), task$replicate)
    truth <- makeTruth(process, scenario$type, scenario$params, "low")
    sim <- simulateHsmm(truth, nCurves[3], initialType = "aggregateStart", seed = seed)
    geoStart <- startModel(sim$x, process, "geo", seed)

    scored <- scoreMethods(list(
      "THMM (funHMM)" = function() { fit <- fitFunhmm(sim$x, geoStart, process); if (!is.null(fit)) fromThmmFit(fit) },
      "THSMM, npg(10)" = function() fromPackageFit(fitThsmm(startModel(sim$x, process, "npg", seed), sim$x), sim$x),
      # The true family of each state, assigned in drift / mean-level order (the
      # order of both the truth and initializeHsmmEM()'s starting states).
      "THSMM, true family" = function() fromPackageFit(fitThsmm(startModel(sim$x, process, scenario$type, seed), sim$x), sim$x)
    ), truth, sim)
    scored[c("rows", "pmf")]
  }

  tasks <- expand.grid(replicate = seq_len(nReps), scenario = names(nonGeometricScenarios), process = names(processes),
                       stringsAsFactors = FALSE)[, 3:1]
  result3 <- runTasks(tasks, section3Replicate, 3)
  rows <- result3$rows
  ari <- rows[rows$parameter == "ARI", c("process", "scenario", "replicate", "method", "estimate")]

  methods <- intersect(names(methodColours), unique(ari$method))
  settings <- unique(ari[, c("process", "scenario")])
  ariTable <- do.call(rbind, lapply(seq_len(nrow(settings)), function(s) {
    g <- merge(settings[s, ], ari)
    cells <- vapply(methods, function(m) meanSd(g$estimate[g$method == m]), character(1))
    # Failed fits are left out of the means and of the paired comparisons below.
    failed <- merge(settings[s, ], rows[rows$parameter == "failed", ])
    cbind(settings[s, ], t(cells), failures = if (nrow(failed)) paste(table(failed$method), names(table(failed$method)), collapse = "; ") else "")
  }))
  rownames(ariTable) <- NULL

  # Paired differences (same data set) of each THSMM against the THMM.
  wide <- reshape(ari, idvar = c("process", "scenario", "replicate"), timevar = "method", direction = "wide")
  names(wide) <- sub("^estimate\\.", "", names(wide))
  comparisons <- expand.grid(thsmm = c("THSMM, npg(10)", "THSMM, true family"),
                             thmm = "THMM (funHMM)", stringsAsFactors = FALSE)
  pairedTable <- do.call(rbind, lapply(seq_len(nrow(settings)), function(s) {
    g <- merge(settings[s, ], wide)
    do.call(rbind, lapply(seq_len(nrow(comparisons)), function(c) {
      difference <- g[[comparisons$thsmm[c]]] - g[[comparisons$thmm[c]]]
      difference <- difference[!is.na(difference)]
      cbind(settings[s, ], data.frame(
        comparison = paste(comparisons$thsmm[c], "-", comparisons$thmm[c]),
        `mean ARI gain` = mean(difference), `se` = sd(difference) / sqrt(length(difference)),
        `share THSMM better` = mean(difference > 0),
        `Wilcoxon p` = if (length(difference) > 1 && any(difference != 0)) suppressWarnings(wilcox.test(difference)$p.value) else NA,
        check.names = FALSE))
    }))
  }))
  rownames(pairedTable) <- NULL

  saveSection(3, result3, list(
    ari = list(table = ariTable, caption = paste(
      "Viterbi ARI under non-geometric dwell times, mean (sd) over replicates. The THMM is fitted",
      "with funHMM.")),
    paired = list(table = pairedTable, caption = paste(
      "Paired ARI differences, THSMM minus THMM, on the same data sets: mean, standard error,",
      "share of replicates where the THSMM decodes better, and Wilcoxon signed-rank p-value."))
  ))

  cat("\n")
  print(ariTable, row.names = FALSE)
  cat("\n")
  print(format(pairedTable, digits = 3), row.names = FALSE)

  pdf(file.path(outputDir, "section3.pdf"), width = 11, height = 4.5)
  par(mfrow = c(1, length(nonGeometricScenarios)), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (name in names(processes)) {
    for (scenarioName in names(nonGeometricScenarios)) {
      g <- ari[ari$process == name & ari$scenario == scenarioName, ]
      boxplot(estimate ~ factor(method, methods), data = g, col = methodColours[methods], names = seq_along(methods),
              xlab = "method", ylab = "Viterbi ARI", main = scenarioName, ylim = c(min(0, g$estimate, na.rm = TRUE), 1))
      if (scenarioName == names(nonGeometricScenarios)[1]) legend("bottomleft", paste(seq_along(methods), methods), bty = "n", cex = 0.7)
    }
    mtext(paste("Section 3, Viterbi ARI:", processes[[name]]$label), outer = TRUE, font = 2)
  }
  invisible(dev.off())
}

cat(sprintf("\nResults written to %s\n", outputDir))
