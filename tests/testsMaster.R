# Master test suite for the THSMM paper: THMM vs THSMM.
#
#   Sections 0-3   nonparametric L2 mean-curve emissions (from tests/l2Tests.R)
#   Section 4      BM / OU / fBM emissions under non-geometric dwell times
#                  (Section 3 of tests/paperSims.R)
#
# Run from the package root:
#   Rscript tests/testsMaster.R            the sections chosen in "Run settings"
#   Rscript tests/testsMaster.R 0 4        selected sections (0: data plots only, seconds)
# The sections, the number of replicates and the cores are set in the "Run
# settings" block below: start with nReps = 1 for a proof of concept, then switch
# to 50. Environment variables override that block when set: TESTS_REPS
# (replicates), TESTS_CORES (parallel workers).
#
# ---------------------------------------------------------------------------
# Sections 0-3: L2 mean curves
# ---------------------------------------------------------------------------
#
# Data: K = 3 states, each emitting curves on 100 equispaced points of [0, 1]:
# the state's mean curve plus smooth Brownian-bridge noise (truncated
# Karhunen-Loeve expansion, 16 terms, coefficient sd 0.4; the noise of the
# paper's Sec. 5.4 and of funHMM::rbridge()). The mean curves are a times
#   sin(2 pi tau),   2 tau - 1,   cos(2 pi tau) + 1,
# with a = 0.15 / 0.2 / 0.3 for low / medium / high separation (nearest-mean
# classification of single curves at the truth is about 0.90 / 0.965 / 0.999
# correct). Embedded jumps are uniform and the chain starts uniformly.
#
# Emission functional: log b_j(O) = -scale / 2 * sum_k (O(tau_k) - h_j(tau_k))^2
# with scale = 1 (the package's "nonpar" family with norm "L2", as in
# viterbiTest.R). funHMM writes the same functional as a Riemann sum divided by
# sigma^2, so it is fitted with sigma = 1 / sqrt(len * scale): the THMM
# (funHMM) and the THSMM with geometric dwell then evaluate identical
# log-likelihoods (checked: equal EM iterates from a common start).
#
# Fits, all from the same k-means start:
#   THMM (funHMM)              the CRAN package funHMM, thmm(type = "nonpar", norm = "L2")
#   THSMM, geometric           geometric dwell (the THMM inside this package)
#   THSMM, negative binomial   1 + NB(size, prob) dwell, which contains the
#                              geometric as size = 1: geo(p) = 1 + NB(1, p)
#   THSMM, shifted Poisson     shift + Poisson(lambda) dwell
#   THSMM, nonparametric       npg dwell: free probabilities for dwell times 1 to
#                              npgSize - 1 plus a geometric tail (no shape assumed)
#   THSMM, true family         the generating dwell family
# Integer shifts are profiled over shiftCandidates by fitHsmmShiftGrid().
#
# Section 0 -- plots of the synthetic data sets (replicate 1 of every setting):
# curves coloured by true state with the true mean curves / mean functions, and
# the hidden state path of the first 60 curves. data.pdf covers Sections 1-3,
# data-section4.pdf covers Section 4.
#
# Section 1 -- recovering a geometric sojourn. Truth: geometric dwell with leave
# probabilities 0.4 / 0.25 / 0.15. Fitted with the THMM, the geometric THSMM and
# the negative binomial THSMM (which should find size near 1). Output: dwell
# distribution plots and a table of estimates.
#
# Section 2 -- recovering non-geometric sojourns. Truth: shifted Poisson, and
# negative binomial dwell times. Every data set is fitted with three THSMM dwell
# families -- shifted Poisson, negative binomial and nonparametric (npg) -- so
# each truth is fitted once by its own family and twice by others (the THMM's
# implied geometric dwell is drawn for contrast). Output: dwell distribution
# plots, a table of estimates, and a table of dwell-recovery error by family.
#
# Section 3 -- Viterbi ARI, THMM vs THSMM, on the Section 2 truths. Output:
# confusion matrices of the Viterbi paths (summed over replicates after matching
# labels to the truth) and a table of ARIs.
#
# Sections 1-3 run at low, medium and high separation.
#
# ---------------------------------------------------------------------------
# Section 4: BM / OU / fBM emissions (paperSims.R Section 3)
# ---------------------------------------------------------------------------
#
# Viterbi ARI, THMM vs THSMM, when dwell times are not geometric, with the
# three emission processes of Kashlak, Loliencar and Heo (2023, Sec. 5) at low
# separation (500 curves on 101 points of [0, 1], O(0) = 0):
#   BM    Brownian motion with drift c = (-2, 0, 2)
#   OU    Ornstein-Uhlenbeck dO = (b0 - b1 O) dtau + dW, mu = (-1, 0, 1), c = 2
#         (b0 = c mu, b1 = c)
#   fBM   fractional BM, hurst 0.75, with linear drift c = (-2, 0, 2)
# Dwell truths: shifted Poisson, peaked negative binomial, Weibull with
# increasing hazard, and a mix of the three families across states. Methods:
#   THMM (funHMM)         funHMM::thmm() with its automatic sigma; it gets the
#                         curves without the O(0) column (its grid convention)
#   THSMM, geometric      the THMM inside this package
#   THSMM, nonparametric  npg(10) dwell
#   THSMM, true family    the generating family per state (shifts profiled)
# (paperSims.R used the archived thmm1.R as its THMM; funHMM is its published
# version.) Seeds are those of paperSims.R Section 3, so the data sets are
# identical. Output: ARI table, paired ARI gains, ARI boxplots (section4.pdf).
#
# Output goes to tests/output/testsMaster/ (tests/output/testsMaster-poc/ when
# nReps = 1).


# ---------------------------------------------------------------------------
# Run settings: edit these
# ---------------------------------------------------------------------------

sections <- 0      # which sections to run (0 = data plots)
nReps <- 50           # replicates per setting: 1 = proof of concept, 50 for the paper
nCores <- parallel::detectCores() - 1
# TRUE redraws the tables and plots from the saved results (sectionN.rds in the
# output folder) without refitting anything; keep nReps at the value of the run
# that produced them, so the right folder is read.
reuseResults <- TRUE 

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

# Sections given on the command line replace the ones in "Run settings".
if (length(commandArgs(trailingOnly = TRUE)) > 0L) sections <- as.integer(commandArgs(trailingOnly = TRUE))

nReps <- as.integer(Sys.getenv("TESTS_REPS", nReps))
nCores <- max(1L, as.integer(Sys.getenv("TESTS_CORES", nCores)))

outputDir <- file.path("tests", "output", if (nReps == 1) "testsMaster-poc" else "testsMaster")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

K <- 3
nPoints <- 100
grid <- seq(0, 1, length.out = nPoints)
emissionScale <- 1
emissionControl <- list(grid = grid, norm = "L2", scale = emissionScale)
funHmmSigma <- 1 / sqrt(nPoints * emissionScale)

shapes <- rbind(sin(2 * pi * grid), 2 * grid - 1, cos(2 * pi * grid) + 1)
separations <- c(low = 0.15, medium = 0.2, high = 0.3)

nCurves <- c(1000, 1000, 1000)      # per section
parametricDwellSize <- 30         # exact dwell support before the geometric tail
npgSize <- 15                     # npg: 14 free head probabilities + geometric tail
shiftCandidates <- 1:5
maxDwell <- 300                   # support for mean dwells and TV distances
maxPlot <- 15                     # dwell support shown in the figures
emTol <- 1e-8
emMaxIter <- 2000

withSize <- function(...) c(..., dwellSize = parametricDwellSize)

geometricLeave <- c(0.4, 0.25, 0.15)                                 # mean dwells 2.5, 4, 6.7
geometricTruth <- list(type = rep("geo", K), params = lapply(geometricLeave, function(p) c(prob = p)))

dwellTruths <- list(
  `Shifted Poisson` = list(type = rep("shiftp", K),
                           params = list(withSize(lambda = 1, shift = 3), withSize(lambda = 2, shift = 2),
                                         withSize(lambda = 1, shift = 4))),               # means 4, 4, 5
  `Negative binomial` = list(type = rep("nb", K),
                             params = list(withSize(size = 3, prob = 0.5), withSize(size = 5, prob = 0.5),
                                           withSize(size = 10, prob = 0.6)))              # 1 + NB: means 4, 6, 7.7
)

# High-contrast colours (black, blue, purple, vermillion, ochre, green) and a
# distinct point shape per method, so lines stay distinguishable in print and
# in greyscale.
methodColours <- c("THMM (funHMM)" = "#000000", "THSMM, geometric" = "#0072B2",
                   "THSMM, shifted Poisson" = "#AA3377", "THSMM, negative binomial" = "#D55E00",
                   "THSMM, nonparametric" = "#997700", "THSMM, true family" = "#009E73")
methodShapes <- c("THMM (funHMM)" = 15, "THSMM, geometric" = 17, "THSMM, shifted Poisson" = 18,
                  "THSMM, negative binomial" = 16, "THSMM, nonparametric" = 25, "THSMM, true family" = 8)
methodOrder <- names(methodColours)


# ---------------------------------------------------------------------------
# Models and fits
# ---------------------------------------------------------------------------

allPermutations <- local({
  build <- function(n) {
    if (n == 1L) return(matrix(1L, 1, 1))
    smaller <- build(n - 1L)
    do.call(rbind, lapply(seq_len(n), function(first) cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))))
  }
  build(K)
})

# match[i] = estimated state that best matches true state i (maximizes matched counts).
matchStates <- function(estimated, truth) {
  tab <- table(factor(estimated, seq_len(K)), factor(truth, seq_len(K)))
  allPermutations[which.max(apply(allPermutations, 1, function(p) sum(tab[cbind(p, seq_len(K))]))), ]
}

makeTruth <- function(dwell, separation) {
  omega <- matrix(1 / (K - 1), K, K)
  diag(omega) <- 0
  buildHsmm(omega, dwell$params, dwell$type, separations[[separation]] * shapes, "nonpar",
            delta = rep(1 / K, K), stationary = FALSE, emissionControl = emissionControl)
}

startModel <- function(x, sojournType, seed) {
  initializeHsmmEM(x, K, sojournType, "nonpar", emissionControl,
                   dwellSize = if (sojournType == "npg") npgSize else parametricDwellSize,
                   stationary = FALSE, seed = seed)$model
}

fitThsmm <- function(start, x) {
  if (any(start$sojournType == "shiftp")) {
    # Integer shifts are held fixed by EM: profile them over a grid.
    shiftValues <- lapply(start$sojournType, function(type) if (type == "shiftp") shiftCandidates)
    return(suppressWarnings(fitHsmmShiftGrid(start, x, shiftValues, stationary = FALSE, fitFunction = fitHsmmEM,
                                             maxCombinations = length(shiftCandidates)^K, tol = emTol, maxIter = emMaxIter)))
  }
  suppressWarnings(fitHsmmEM(start, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
}

thmmMatrix <- function(model) {
  leave <- vapply(model$sojournParameters, function(s) s[["prob"]], numeric(1))
  diag(1 - leave) + leave * model$tpm
}

# Fitted methods in a common form: Viterbi path, mean curves, and per-(estimated)
# state dwell pmf / family / parameters.
fromPackageFit <- function(fit, x) {
  model <- fit$fittedModel
  list(states = hsmmViterbiEM(model, x)$states, means = model$emissionParameters,
       dwellPmf = function(j) getModelDwellPmf(model, j, maxDwell),
       dwellType = function(j) model$sojournType[j],
       dwellParams = function(j) model$sojournParameters[[j]])
}

fitFunHmm <- function(x, start) {
  fit <- suppressWarnings(funHMM::thmm(x, K, type = "nonpar", norm = "L2", sigma = funHmmSigma,
                                       par = start$emissionParameters, trans = thmmMatrix(start),
                                       init.prob = start$delta, tol = emTol, max.iter = emMaxIter))
  list(states = fit$states, means = fit$par,
       dwellPmf = function(j) dgeom(seq_len(maxDwell) - 1, 1 - fit$trans[j, j]),
       dwellType = function(j) "geo",
       dwellParams = function(j) c(prob = 1 - fit$trans[j, j]))
}

runMethods <- function(x, seed, which, trueType) {
  geoStart <- startModel(x, "geo", seed)
  methods <- list(
    "THMM (funHMM)" = function() fitFunHmm(x, geoStart),
    "THSMM, geometric" = function() fromPackageFit(fitThsmm(geoStart, x), x),
    "THSMM, shifted Poisson" = function() fromPackageFit(fitThsmm(startModel(x, "shiftp", seed), x), x),
    "THSMM, negative binomial" = function() fromPackageFit(fitThsmm(startModel(x, "nb", seed), x), x),
    "THSMM, nonparametric" = function() fromPackageFit(fitThsmm(startModel(x, "npg", seed), x), x),
    "THSMM, true family" = function() fromPackageFit(fitThsmm(startModel(x, trueType, seed), x), x)
  )
  lapply(methods[which], function(f) tryCatch(f(), error = function(e) list(error = conditionMessage(e))))
}


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

resultRow <- function(state, parameter, truth, estimate) {
  data.frame(state = state, parameter = parameter, truth = truth, estimate = estimate)
}

# Dwell parameters of the truth in a fitted family, when the family contains it:
# geo(p) is 1 + NB(size = 1, prob = p).
trueDwellIn <- function(fittedType, truth, i) {
  trueType <- truth$sojournType[i]
  params <- truth$sojournParameters[[i]]
  if (fittedType == trueType) return(params[setdiff(names(params), "dwellSize")])
  if (fittedType == "nb" && trueType == "geo") return(c(size = 1, prob = params[["prob"]]))
  NULL
}

scoreMethod <- function(method, fitted, truth, trueStates) {
  match <- matchStates(fitted$states, trueStates)
  matched <- order(match)[fitted$states]   # decoded labels renamed to their matched true state
  rows <- list(resultRow(NA, "ARI", 1, funHMM::ari(fitted$states, trueStates)))
  dwellTimes <- seq_len(maxDwell)
  pmfRows <- list()

  for (i in seq_len(K)) {
    j <- match[i]
    rows[[length(rows) + 1L]] <- resultRow(i, "mean curve RMSE", 0,
                                           sqrt(mean((fitted$means[j, ] - truth$emissionParameters[i, ])^2)))
    truePmf <- getModelDwellPmf(truth, i, maxDwell)
    estPmf <- fitted$dwellPmf(j)
    rows[[length(rows) + 1L]] <- resultRow(i, "mean dwell", sum(dwellTimes * truePmf), sum(dwellTimes * estPmf))
    rows[[length(rows) + 1L]] <- resultRow(i, "TV dwell", 0, 0.5 * sum(abs(truePmf - estPmf)))
    trueDwell <- trueDwellIn(fitted$dwellType(j), truth, i)
    for (p in names(trueDwell)) {
      rows[[length(rows) + 1L]] <- resultRow(i, p, trueDwell[[p]], fitted$dwellParams(j)[[p]])
    }
    pmfRows[[i]] <- data.frame(method = method, state = i, time = seq_len(maxPlot), prob = estPmf[seq_len(maxPlot)])
  }

  confusion <- as.data.frame(table(truth = factor(trueStates, seq_len(K)), viterbi = factor(matched, seq_len(K))),
                             responseName = "count")
  list(rows = cbind(method = method, do.call(rbind, rows)), pmf = do.call(rbind, pmfRows),
       confusion = cbind(method = method, confusion))
}

# Every data set has its own seed, so any replicate can be rebuilt exactly.
datasetSeed <- function(section, dwellName, separation, replicate) {
  as.integer(1e6 * section + 1e4 * match(separation, names(separations)) +
               1e2 * match(dwellName, c("Geometric", names(dwellTruths))) + replicate)
}

replicateFunction <- function(section, which) {
  function(task) {
    dwell <- if (task$dwell == "Geometric") geometricTruth else dwellTruths[[task$dwell]]
    seed <- datasetSeed(section, task$dwell, task$separation, task$replicate)
    truth <- makeTruth(dwell, task$separation)
    sim <- simulateHsmm(truth, nCurves[section], initialType = "aggregateStart", seed = seed)
    fits <- runMethods(sim$x, seed, which, dwell$type[1])
    scored <- lapply(names(fits), function(m) {
      if (!is.null(fits[[m]]$error)) {
        return(list(rows = data.frame(method = m, state = NA, parameter = "failed", truth = NA, estimate = 1)))
      }
      scoreMethod(m, fits[[m]], truth, sim$states)
    })
    lapply(c(rows = "rows", pmf = "pmf", confusion = "confusion"), function(field) do.call(rbind, lapply(scored, `[[`, field)))
  }
}

runTasks <- function(truths, section, which) {
  tasks <- expand.grid(replicate = seq_len(nReps), separation = names(separations), dwell = truths,
                       stringsAsFactors = FALSE)[, 3:1]
  started <- Sys.time()
  f <- replicateFunction(section, which)
  outputs <- mclapply(seq_len(nrow(tasks)), function(t) {
    tryCatch(f(tasks[t, , drop = FALSE]), error = function(e) list(error = conditionMessage(e)))
  }, mc.cores = nCores, mc.preschedule = FALSE)

  failed <- vapply(outputs, function(o) !is.list(o) || !is.null(o$error), logical(1))
  if (any(failed)) cat(sprintf("  %d of %d replicates failed, e.g.: %s\n", sum(failed), length(outputs),
                               as.character(outputs[[which(failed)[1]]]$error %||% outputs[[which(failed)[1]]])))
  attach <- function(field) do.call(rbind, lapply(which(!failed), function(t) {
    if (is.null(outputs[[t]][[field]])) return(NULL)
    cbind(tasks[t, , drop = FALSE], outputs[[t]][[field]], row.names = NULL)
  }))
  cat(sprintf("  Section %d: %d replicates in %.1f minutes on %d cores\n", section, nrow(tasks),
              as.numeric(difftime(Sys.time(), started, units = "mins")), nCores))
  result <- list(rows = attach("rows"), pmf = attach("pmf"), confusion = attach("confusion"))
  failures <- result$rows[result$rows$parameter == "failed", ]
  if (nrow(failures)) cat(sprintf("  %d failed fits: %s\n", nrow(failures), paste(unique(failures$method), collapse = ", ")))
  result
}


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

separationOrder <- function(s) match(s, names(separations))

# Bias / sd / RMSE of every (dwell truth, separation, method, state, parameter).
summarizeEstimates <- function(rows) {
  rows <- rows[!is.na(rows$truth) & rows$parameter != "ARI", ]
  by <- c("dwell", "separation", "method", "state", "parameter")
  groups <- split(rows, rows[, by], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    error <- g$estimate - g$truth
    cbind(g[1, by, drop = FALSE],
          data.frame(true = g$truth[1], mean = mean(g$estimate, na.rm = TRUE), bias = mean(error, na.rm = TRUE),
                     sd = sd(g$estimate, na.rm = TRUE), rmse = sqrt(mean(error^2, na.rm = TRUE)),
                     n = sum(!is.na(g$estimate))))
  }))
  rownames(out) <- NULL
  out[order(out$dwell, separationOrder(out$separation), match(out$method, methodOrder), out$state), ]
}

writeLatex <- function(table, file, caption, digits = 3) {
  escape <- function(s) gsub("([&%#_$])", "\\\\\\1", s)
  cells <- vapply(table, function(column) {
    if (is.numeric(column)) formatC(column, digits = digits, format = "g") else escape(as.character(column))
  }, character(nrow(table)))
  cells <- matrix(cells, nrow = nrow(table))
  align <- paste(ifelse(vapply(table, is.numeric, logical(1)), "r", "l"), collapse = "")
  writeLines(c("\\begin{table}[ht]", "\\centering", sprintf("\\caption{%s}", caption),
               sprintf("\\begin{tabular}{%s}", align), "\\hline",
               paste(paste(escape(names(table)), collapse = " & "), "\\\\"), "\\hline",
               paste(apply(cells, 1, paste, collapse = " & "), "\\\\"),
               "\\hline", "\\end{tabular}", "\\end{table}"), file)
}

writeTable <- function(table, name, caption) {
  write.csv(table, file.path(outputDir, paste0(name, ".csv")), row.names = FALSE)
  writeLatex(table, file.path(outputDir, paste0(name, ".tex")), caption)
}

# One page per (truth, separation): K panels of dwell pmfs, truth (grey bars),
# fitted median per method (lines) and the 10-90% band of `bandMethod` (one
# method, or a vector of methods named by dwell truth).
dwellPlots <- function(result, truths, bandMethods, title) {
  for (truthName in truths) {
    bandMethod <- if (length(bandMethods) > 1) bandMethods[[truthName]] else bandMethods
    dwell <- if (truthName == "Geometric") geometricTruth else dwellTruths[[truthName]]
    for (separation in names(separations)) {
      truth <- makeTruth(dwell, separation)
      p <- result$pmf[result$pmf$dwell == truthName & result$pmf$separation == separation, ]
      methods <- intersect(methodOrder, unique(p$method))
      par(mfrow = c(1, K), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0), las = 1)
      for (i in seq_len(K)) {
        truePmf <- getModelDwellPmf(truth, i, maxPlot)
        quantiles <- lapply(methods, function(m) {
          d <- p[p$method == m & p$state == i, ]
          sapply(split(d$prob, d$time), quantile, c(0.1, 0.5, 0.9), na.rm = TRUE)
        })
        names(quantiles) <- methods
        upper <- max(truePmf, unlist(lapply(quantiles, function(q) q[3, ])))
        shown <- setdiff(names(dwell$params[[i]]), "dwellSize")
        plot(seq_len(maxPlot), truePmf, type = "h", lwd = 7, col = "grey80", ylim = c(0, upper * 1.05),
             xlab = "dwell time (curves)", ylab = "probability", font.main = 1, cex.main = 0.95,
             main = sprintf("State %d: %s", i, paste(shown, signif(dwell$params[[i]][shown], 3), sep = " = ", collapse = ", ")))
        if (bandMethod %in% methods) {
          polygon(c(seq_len(maxPlot), rev(seq_len(maxPlot))), c(quantiles[[bandMethod]][1, ], rev(quantiles[[bandMethod]][3, ])),
                  col = adjustcolor(methodColours[[bandMethod]], 0.2), border = NA)
        }
        for (m in methods) {
          lines(seq_len(maxPlot), quantiles[[m]][2, ], type = "b", lwd = 1.5, pch = methodShapes[[m]], cex = 0.9,
                col = methodColours[[m]], bg = methodColours[[m]])
        }
        if (i == 1) legend("topright", c("truth", methods), col = c("grey80", methodColours[methods]),
                           pt.bg = c(NA, methodColours[methods]), lwd = c(7, rep(1.5, length(methods))),
                           pch = c(NA, methodShapes[methods]), bty = "n", cex = 0.75)
      }
      mtext(sprintf("%s: %s dwell times, %s separation (a = %g)%s", title, tolower(truthName), separation,
                    separations[[separation]], if (nReps > 1) sprintf(", median and 10-90%% band over %d replicates", nReps) else ""),
            outer = TRUE, font = 2)
    }
  }
}


# ---------------------------------------------------------------------------
# Section 4 settings and helpers: BM / OU / fBM emissions
# ---------------------------------------------------------------------------

s4Grid <- seq(0, 1, length.out = 101)
s4Curves <- 500
s4NpgSize <- 10                   # npg: 9 free head probabilities + geometric tail
s4ShiftCandidates <- 1:4

s4OuParams <- function(mu, rate) cbind(b0 = rate * mu, b1 = rate)

# Low separation: single-curve classification at the true parameters is about
# 0.78 correct for each process. States are ordered by drift / mean level, as
# initializeHsmmEM() orders its starting states.
s4Processes <- list(
  BM = list(label = "Brownian motion with drift", emissionType = "bmwd", hurst = 0.5,
            control = list(grid = s4Grid), params = cbind(drift = c(-2, 0, 2))),
  OU = list(label = "Ornstein-Uhlenbeck", emissionType = "ou", hurst = 0.5,
            control = list(grid = s4Grid), params = s4OuParams(c(-1, 0, 1), 2)),
  fBM = list(label = "Fractional BM with drift (H = 0.75)", emissionType = "bmwd", hurst = 0.75,
             control = list(grid = s4Grid, hurst = 0.75), params = cbind(drift = c(-2, 0, 2)))
)

# Non-geometric dwell times with a typical length (mode > 1).
s4Scenarios <- list(
  `Shifted Poisson` = list(type = rep("shiftp", K),
                           params = list(withSize(lambda = 1, shift = 3), withSize(lambda = 2, shift = 3),
                                         withSize(lambda = 1, shift = 4))),                    # means 4, 5, 5
  `Peaked negative binomial` = list(type = rep("nb", K),
                                    params = list(withSize(size = 20, prob = 20 / 23), withSize(size = 20, prob = 20 / 24),
                                                  withSize(size = 20, prob = 20 / 25))),       # means 4, 5, 6
  `Weibull (increasing hazard)` = list(type = rep("weibull", K),
                                       params = list(withSize(shape = 3, scale = 4), withSize(shape = 3, scale = 5),
                                                     withSize(shape = 3, scale = 6))),          # means 4.1, 5, 5.9
  # The true family of each state is assigned in drift / mean-level order (the
  # order of both the truth and initializeHsmmEM()'s starting states).
  `Mixed families` = list(type = c("shiftp", "nb", "weibull"),
                          params = list(withSize(lambda = 1, shift = 3), withSize(size = 20, prob = 20 / 24),
                                        withSize(shape = 3, scale = 6)))
)

# Section 4 is Section 3 of the write-up, so its plot titles say so.
s4PlotLabel <- "Section 3"

s4Methods <- c("THMM (funHMM)", "THSMM, geometric", "THSMM, nonparametric", "THSMM, true family")

s4MakeTruth <- function(process, scenario) {
  omega <- matrix(1 / (K - 1), K, K)
  diag(omega) <- 0
  buildHsmm(omega, scenario$params, scenario$type, process$params, process$emissionType,
            delta = rep(1 / K, K), stationary = FALSE, emissionControl = process$control)
}

# The seeds of paperSims.R Section 3 (section 3, process index, scenario index,
# replicate), so the data sets are identical to that study's.
s4Seed <- function(processName, scenarioName, replicate) {
  as.integer(1e6 * 3 + 1e4 * match(processName, names(s4Processes)) +
               1e2 * match(scenarioName, names(s4Scenarios)) + replicate)
}

s4Simulate <- function(processName, scenarioName, replicate) {
  truth <- s4MakeTruth(s4Processes[[processName]], s4Scenarios[[scenarioName]])
  sim <- simulateHsmm(truth, s4Curves, initialType = "aggregateStart",
                      seed = s4Seed(processName, scenarioName, replicate))
  c(sim, list(truth = truth))
}

s4StartModel <- function(x, process, sojournType, seed) {
  sojournType <- rep_len(sojournType, K)
  initializeHsmmEM(x, K, sojournType, process$emissionType, process$control,
                   dwellSize = ifelse(sojournType == "npg", s4NpgSize, parametricDwellSize),
                   stationary = FALSE, seed = seed)$model
}

s4FitThsmm <- function(start, x) {
  if (any(start$sojournType == "shiftp")) {
    shiftValues <- lapply(start$sojournType, function(type) if (type == "shiftp") s4ShiftCandidates)
    return(suppressWarnings(fitHsmmShiftGrid(start, x, shiftValues, stationary = FALSE, fitFunction = fitHsmmEM,
                                             maxCombinations = length(s4ShiftCandidates)^K, tol = emTol,
                                             maxIter = emMaxIter)))
  }
  suppressWarnings(fitHsmmEM(start, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
}

# funHMM's THMM from the package's geometric start: drift vector, or (mean,
# rate) for OU; curves without the O(0) column (funHMM's grid is k / len).
s4FitFunHmm <- function(x, start, process) {
  params <- start$emissionParameters
  par <- if (process$emissionType == "ou") {
    rate <- pmax(params[, 2], 1e-3)   # a start at the rate = 0 boundary has no mean
    cbind(mean = params[, 1] / rate, rate = rate)
  } else {
    params[, 1]
  }
  fit <- suppressWarnings(funHMM::thmm(x[, -1], K, type = process$emissionType, hurst = process$hurst, par = par,
                                       trans = thmmMatrix(start), init.prob = start$delta, tol = emTol,
                                       max.iter = emMaxIter))
  fit$states
}

s4Replicate <- function(task) {
  process <- s4Processes[[task$process]]
  scenario <- s4Scenarios[[task$scenario]]
  seed <- s4Seed(task$process, task$scenario, task$replicate)
  sim <- s4Simulate(task$process, task$scenario, task$replicate)
  x <- sim$x
  geoStart <- s4StartModel(x, process, "geo", seed)
  decode <- function(fit) hsmmViterbiEM(fit$fittedModel, x)$states
  methods <- list(
    "THMM (funHMM)" = function() s4FitFunHmm(x, geoStart, process),
    "THSMM, geometric" = function() decode(s4FitThsmm(geoStart, x)),
    "THSMM, nonparametric" = function() decode(s4FitThsmm(s4StartModel(x, process, "npg", seed), x)),
    "THSMM, true family" = function() decode(s4FitThsmm(s4StartModel(x, process, scenario$type, seed), x))
  )
  do.call(rbind, lapply(names(methods), function(m) {
    states <- tryCatch(methods[[m]](), error = function(e) NULL)
    if (is.null(states)) return(data.frame(method = m, parameter = "failed", estimate = 1))
    data.frame(method = m, parameter = "ARI", estimate = funHMM::ari(states, sim$states))
  }))
}

s4RunTasks <- function() {
  tasks <- expand.grid(replicate = seq_len(nReps), scenario = names(s4Scenarios), process = names(s4Processes),
                       stringsAsFactors = FALSE)[, 3:1]
  started <- Sys.time()
  outputs <- mclapply(seq_len(nrow(tasks)), function(t) {
    tryCatch(s4Replicate(tasks[t, , drop = FALSE]), error = function(e) conditionMessage(e))
  }, mc.cores = nCores, mc.preschedule = FALSE)
  failed <- !vapply(outputs, is.data.frame, logical(1))
  if (any(failed)) cat(sprintf("  %d of %d replicates failed, e.g.: %s\n", sum(failed), length(outputs),
                               as.character(outputs[[which(failed)[1]]])))
  cat(sprintf("  Section 4: %d replicates in %.1f minutes on %d cores\n", nrow(tasks),
              as.numeric(difftime(Sys.time(), started, units = "mins")), nCores))
  do.call(rbind, lapply(which(!failed), function(t) cbind(tasks[t, , drop = FALSE], outputs[[t]], row.names = NULL)))
}

# Mean function E O(tau) of each state: c tau for BM and fBM, and
# mu (1 - exp(-c tau)) for OU started at 0.
s4MeanFunctions <- function(process) {
  params <- process$params
  if (process$emissionType == "ou") return(t(sapply(seq_len(K), function(i) params[i, 1] / params[i, 2] * (1 - exp(-params[i, 2] * s4Grid)))))
  t(sapply(seq_len(K), function(i) params[i, 1] * s4Grid))
}

s4ParameterText <- function(process) {
  params <- process$params
  if (process$emissionType == "ou") {
    return(sprintf("mu = (%s), c = %g", paste(params[, 1] / params[, 2], collapse = ", "), params[1, 2]))
  }
  sprintf("c = (%s)", paste(params[, 1], collapse = ", "))
}

# One page per data set: curves coloured by true state with the state mean
# functions, and the hidden state path of the first `nShown` curves.
s4PlotData <- function(sim, process, main, nShown = 60) {
  stateColours <- c("#0072B2", "#D55E00", "#009E73")   # Okabe-Ito
  layout(matrix(1:2, 2), heights = c(3, 1))
  par(mar = c(4, 4, 2.5, 1), las = 1)
  matplot(s4Grid, t(sim$x), type = "l", lty = 1, lwd = 0.5, col = adjustcolor(stateColours[sim$states], 0.3),
          xlab = expression(tau), ylab = expression(O(tau)), main = main, font.main = 1)
  matlines(s4Grid, t(s4MeanFunctions(process)), lty = 1, lwd = 2.5, col = stateColours)
  legend("topleft", c(paste("State", seq_len(K)), s4ParameterText(process)), col = c(stateColours, NA),
         lwd = c(rep(2.5, K), NA), bty = "n", cex = 0.8)
  par(mar = c(4, 4, 0.5, 1))
  plot(seq_len(nShown), sim$states[seq_len(nShown)], type = "s", yaxt = "n", ylim = c(0.5, K + 0.5),
       xlab = "time (curve index)", ylab = "state")
  axis(2, at = seq_len(K))
  layout(1)
}

cat(sprintf("testsMaster: sections %s, %d replicates, %d cores\n", paste(sections, collapse = ", "), nReps, nCores))


# ===========================================================================
# Section 0: plots of the synthetic data sets
# ===========================================================================

if (0 %in% sections) {
  cat("\nSection 0: data plots\n")
  stateColours <- c("#0072B2", "#D55E00", "#009E73")   # Okabe-Ito
  nShown <- 60

  pdf(file.path(outputDir, "data.pdf"), width = 8, height = 6)
  for (section in 1:3) {
    for (dwellName in if (section == 1) "Geometric" else names(dwellTruths)) {
      dwell <- if (dwellName == "Geometric") geometricTruth else dwellTruths[[dwellName]]
      for (separation in names(separations)) {
        # Replicate 1, exactly as the fits in that section see it.
        truth <- makeTruth(dwell, separation)
        sim <- simulateHsmm(truth, nCurves[section], initialType = "aggregateStart",
                            seed = datasetSeed(section, dwellName, separation, 1))

        layout(matrix(1:2, 2), heights = c(3, 1))
        par(mar = c(4, 4, 2.5, 1), las = 1)
        matplot(grid, t(sim$x), type = "l", lty = 1, lwd = 0.5, col = adjustcolor(stateColours[sim$states], 0.35),
                xlab = expression(tau), ylab = expression(O(tau)), font.main = 1,
                main = sprintf("Section %d: %s dwell, %s separation (a = %g), %d curves", section, tolower(dwellName),
                               separation, separations[[separation]], nCurves[section]))
        matlines(grid, t(truth$emissionParameters), lty = 1, lwd = 2.5, col = stateColours)
        legend("topright", paste("State", seq_len(K)), col = stateColours, lwd = 2.5, bty = "n", cex = 0.8)
        par(mar = c(4, 4, 0.5, 1))
        plot(seq_len(nShown), sim$states[seq_len(nShown)], type = "s", yaxt = "n", ylim = c(0.5, K + 0.5),
             xlab = "time (curve index)", ylab = "state")
        axis(2, at = seq_len(K))
        layout(1)
      }
    }
  }
  invisible(dev.off())
  cat(sprintf("  Wrote %s\n", file.path(outputDir, "data.pdf")))

  # Section 4 data sets (replicate 1 of every process x dwell scenario).
  pdf(file.path(outputDir, "data-section4.pdf"), width = 8, height = 6)
  for (processName in names(s4Processes)) {
    for (scenarioName in names(s4Scenarios)) {
      sim <- s4Simulate(processName, scenarioName, 1)
      s4PlotData(sim, s4Processes[[processName]], sprintf("%s: %s, %s dwell, %d curves", s4PlotLabel,
                                                          s4Processes[[processName]]$label, tolower(scenarioName), s4Curves))
    }
  }
  invisible(dev.off())
  cat(sprintf("  Wrote %s\n", file.path(outputDir, "data-section4.pdf")))
}


# ===========================================================================
# Section 1: recovering a geometric sojourn distribution
# ===========================================================================

if (1 %in% sections) {
  cat("\nSection 1: geometric sojourn recovery\n")
  if (reuseResults) {
    result1 <- readRDS(file.path(outputDir, "section1.rds"))
  } else {
    result1 <- runTasks("Geometric", 1, c("THMM (funHMM)", "THSMM, geometric", "THSMM, negative binomial"))
    saveRDS(result1, file.path(outputDir, "section1.rds"))
  }
  write.csv(result1$rows, file.path(outputDir, "section1-replicates.csv"), row.names = FALSE)

  estimates <- summarizeEstimates(result1$rows)
  writeTable(estimates, "section1-estimates", paste(
    "Geometric truth: estimates of matched states (leave probability prob, negative binomial size and prob,",
    "mean dwell, TV distance between fitted and true dwell pmfs, RMSE of the fitted mean curve),",
    "by separation; geo(p) = 1 + NB(size = 1, prob = p)."))
  cat("\n")
  print(format(estimates[estimates$parameter %in% c("prob", "size", "mean dwell", "TV dwell"), -1], digits = 3), row.names = FALSE)

  pdf(file.path(outputDir, "section1-dwell.pdf"), width = 10, height = 3.8)
  dwellPlots(result1, "Geometric", "THSMM, negative binomial", "Section 1")
  invisible(dev.off())
}


# ===========================================================================
# Section 2: recovering shifted Poisson and negative binomial sojourns
# ===========================================================================

if (2 %in% sections) {
  cat("\nSection 2: non-geometric sojourn recovery\n")
  fittedFamilies <- c("THSMM, shifted Poisson", "THSMM, negative binomial", "THSMM, nonparametric")
  if (reuseResults) {
    result2 <- readRDS(file.path(outputDir, "section2.rds"))
  } else {
    result2 <- runTasks(names(dwellTruths), 2, c("THMM (funHMM)", fittedFamilies))
    saveRDS(result2, file.path(outputDir, "section2.rds"))
  }
  write.csv(result2$rows, file.path(outputDir, "section2-replicates.csv"), row.names = FALSE)

  estimates <- summarizeEstimates(result2$rows)
  writeTable(estimates, "section2-estimates", paste(
    "Non-geometric truths: dwell and mean-curve estimates of matched states, by separation."))
  cat("\n")
  print(format(estimates[estimates$parameter %in% c("lambda", "shift", "size", "prob"), ], digits = 3), row.names = FALSE)

  # Dwell-recovery error by fitted family: mean TV distance (over states and
  # replicates) and RMSE of the state mean dwells.
  rows <- result2$rows
  recovery <- do.call(rbind, lapply(split(rows, rows[, c("dwell", "separation", "method")], drop = TRUE), function(g) {
    tv <- g$estimate[g$parameter == "TV dwell"]
    meanDwell <- g[g$parameter == "mean dwell", ]
    data.frame(dwell = g$dwell[1], separation = g$separation[1], method = g$method[1],
               `mean TV dwell` = mean(tv), `mean dwell RMSE` = sqrt(mean((meanDwell$estimate - meanDwell$truth)^2)),
               ARI = mean(g$estimate[g$parameter == "ARI"]), check.names = FALSE)
  }))
  recovery <- recovery[order(recovery$dwell, separationOrder(recovery$separation), match(recovery$method, methodOrder)), ]
  rownames(recovery) <- NULL
  writeTable(recovery, "section2-recovery", paste(
    "Dwell-recovery error by fitted family: mean TV distance between fitted and true dwell pmfs, RMSE of",
    "the state mean dwells, and Viterbi ARI, by dwell truth and separation."))
  cat("\n")
  print(format(recovery, digits = 3), row.names = FALSE)

  pdf(file.path(outputDir, "section2-dwell.pdf"), width = 10, height = 3.8)
  dwellPlots(result2, names(dwellTruths),
             c(`Shifted Poisson` = "THSMM, shifted Poisson", `Negative binomial` = "THSMM, negative binomial"), "Section 2")
  invisible(dev.off())
}


# ===========================================================================
# Section 3: Viterbi ARI, THMM vs THSMM
# ===========================================================================

if (3 %in% sections) {
  cat("\nSection 3: Viterbi ARI, THMM vs THSMM\n")
  if (reuseResults) {
    result3 <- readRDS(file.path(outputDir, "section3.rds"))
  } else {
    result3 <- runTasks(names(dwellTruths), 3, c("THMM (funHMM)", "THSMM, geometric", "THSMM, true family"))
    saveRDS(result3, file.path(outputDir, "section3.rds"))
  }
  write.csv(result3$rows, file.path(outputDir, "section3-replicates.csv"), row.names = FALSE)

  # ARI table: mean (sd) per truth x separation x method.
  ari <- result3$rows[result3$rows$parameter == "ARI", ]
  cell <- function(v) if (length(v) > 1) sprintf("%.3f (%.3f)", mean(v), sd(v)) else sprintf("%.3f", v)
  ariTable <- do.call(rbind, lapply(split(ari, ari[, c("dwell", "separation")], drop = TRUE), function(g) {
    methods <- intersect(methodOrder, unique(g$method))
    cbind(g[1, c("dwell", "separation")], as.data.frame(t(vapply(methods, function(m) cell(g$estimate[g$method == m]), "")),
                                                        check.names = FALSE))
  }))
  ariTable <- ariTable[order(ariTable$dwell, separationOrder(ariTable$separation)), ]
  rownames(ariTable) <- NULL
  writeTable(ariTable, "section3-ari", paste(
    "Viterbi ARI against the true states", if (nReps > 1) sprintf("(mean (sd) over %d replicates)", nReps) else "(one data set)",
    "by dwell truth and separation. THSMM, geometric is the THMM fitted by the THSMM code."))
  cat("\n")
  print(ariTable, row.names = FALSE)

  # Confusion matrices: true state (rows) vs Viterbi state matched to the truth
  # (columns), summed over replicates.
  confusion <- aggregate(count ~ dwell + separation + method + truth + viterbi, result3$confusion, sum)
  names(confusion)[names(confusion) == "truth"] <- "trueState"
  write.csv(confusion, file.path(outputDir, "section3-confusion.csv"), row.names = FALSE)
  sink(file.path(outputDir, "section3-confusion.txt"))
  cat(sprintf("Viterbi confusion matrices: true state (rows) vs decoded state matched to the truth (columns)%s\n",
              if (nReps > 1) sprintf(", summed over %d replicates", nReps) else ""))
  for (truthName in names(dwellTruths)) for (separation in names(separations)) {
    for (method in intersect(methodOrder, unique(confusion$method))) {
      g <- confusion[confusion$dwell == truthName & confusion$separation == separation & confusion$method == method, ]
      if (nrow(g) == 0L) next
      cat(sprintf("\n%s dwell, %s separation | %s\n", truthName, separation, method))
      print(xtabs(count ~ trueState + viterbi, g))
    }
  }
  sink()
  cat(sprintf("\nConfusion matrices written to %s\n", file.path(outputDir, "section3-confusion.txt")))
}

# ===========================================================================
# Section 4: Viterbi ARI, THMM vs THSMM, BM / OU / fBM emissions
# ===========================================================================

if (4 %in% sections) {
  cat("\nSection 4: Viterbi ARI under non-geometric dwell times, BM / OU / fBM\n")
  if (reuseResults) {
    rows4 <- readRDS(file.path(outputDir, "section4.rds"))
  } else {
    rows4 <- s4RunTasks()
    saveRDS(rows4, file.path(outputDir, "section4.rds"))
  }
  write.csv(rows4, file.path(outputDir, "section4-replicates.csv"), row.names = FALSE)

  ari4 <- rows4[rows4$parameter == "ARI", ]
  cell <- function(v) if (length(v) > 1) sprintf("%.3f (%.3f)", mean(v), sd(v)) else if (length(v)) sprintf("%.3f", v) else NA
  settings4 <- unique(rows4[, c("process", "scenario")])
  ariTable4 <- do.call(rbind, lapply(seq_len(nrow(settings4)), function(s) {
    g <- merge(settings4[s, ], ari4)
    failed <- merge(settings4[s, ], rows4[rows4$parameter == "failed", ])
    cbind(settings4[s, ], as.data.frame(t(vapply(s4Methods, function(m) cell(g$estimate[g$method == m]), "")),
                                        check.names = FALSE),
          failures = if (nrow(failed)) paste(table(failed$method), names(table(failed$method)), collapse = "; ") else "")
  }))
  rownames(ariTable4) <- NULL
  writeTable(ariTable4, "section4-ari", paste(
    "Viterbi ARI under non-geometric dwell times, BM / OU / fBM emissions at low separation",
    if (nReps > 1) sprintf("(mean (sd) over %d replicates).", nReps) else "(one data set).",
    "THSMM, geometric is the THMM fitted by the THSMM code; failed fits are excluded."))
  cat("\n")
  print(ariTable4, row.names = FALSE)

  # Paired differences (same data set) of each THSMM against each THMM fit.
  wide4 <- reshape(ari4[, c("process", "scenario", "replicate", "method", "estimate")],
                   idvar = c("process", "scenario", "replicate"), timevar = "method", direction = "wide")
  names(wide4) <- sub("^estimate\\.", "", names(wide4))
  paired4 <- do.call(rbind, lapply(seq_len(nrow(settings4)), function(s) {
    g <- merge(settings4[s, ], wide4)
    do.call(rbind, lapply(c("THSMM, nonparametric", "THSMM, true family"), function(a) {
      do.call(rbind, lapply(c("THMM (funHMM)", "THSMM, geometric"), function(b) {
        difference <- g[[a]] - g[[b]]
        difference <- difference[!is.na(difference)]
        cbind(settings4[s, ], data.frame(
          comparison = paste(a, "-", sub("THSMM, geometric", "THMM (package)", b)),
          `mean ARI gain` = mean(difference),
          se = if (length(difference) > 1) sd(difference) / sqrt(length(difference)) else NA,
          `share THSMM better` = mean(difference > 0),
          `Wilcoxon p` = if (length(difference) > 1 && any(difference != 0)) suppressWarnings(wilcox.test(difference)$p.value) else NA,
          check.names = FALSE))
      }))
    }))
  }))
  rownames(paired4) <- NULL
  writeTable(paired4, "section4-paired", paste(
    "Paired ARI differences, THSMM minus THMM, on the same data sets: mean, standard error, share of",
    "replicates where the THSMM decodes better, and Wilcoxon signed-rank p-value."))
  cat("\n")
  print(format(paired4, digits = 3), row.names = FALSE)

  # ARI boxplots: one page per process, one panel per dwell scenario.
  pdf(file.path(outputDir, "section4.pdf"), width = 11, height = 4.5)
  colours <- methodColours[s4Methods]
  for (processName in names(s4Processes)) {
    par(mfrow = c(1, length(s4Scenarios)), mar = c(3, 4, 2.5, 0.5), oma = c(2.5, 0, 2, 0), las = 1)
    for (scenarioName in names(s4Scenarios)) {
      g <- ari4[ari4$process == processName & ari4$scenario == scenarioName, ]
      boxplot(estimate ~ factor(method, s4Methods), data = g, col = adjustcolor(colours, 0.6), names = seq_along(s4Methods),
              xlab = "", ylab = "Viterbi ARI", main = scenarioName, font.main = 1, ylim = c(min(0, g$estimate), 1))
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", paste(seq_along(s4Methods), s4Methods), fill = adjustcolor(colours, 0.6), horiz = TRUE, bty = "n", cex = 0.8)
    mtext(paste0(s4PlotLabel, ", Viterbi ARI: ", s4Processes[[processName]]$label), side = 3, line = -1.5, font = 2)
  }
  invisible(dev.off())
}

cat(sprintf("\nResults written to %s\n", outputDir))
