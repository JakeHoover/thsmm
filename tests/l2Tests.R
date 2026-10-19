# THMM vs THSMM with nonparametric L2 mean-curve emissions.
#
# Run from the package root:
#   Rscript tests/l2Tests.R            all sections (0-3)
#   Rscript tests/l2Tests.R 1 3        selected sections (0: data plots only, seconds)
# The number of replicates is set in the "Run settings" block below: start with
# nReps = 1 for a proof of concept, then switch to 50. Environment variables
# override that block when set: L2_TESTS_REPS (replicates), L2_TESTS_CORES
# (parallel workers).
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
# Section 0 -- plots of the synthetic data sets (replicate 1 of every section,
# dwell truth and separation): curves coloured by true state with the true mean
# curves, and the hidden state path of the first 60 curves (data.pdf).
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
# Every section runs at low, medium and high separation. Output goes to
# tests/output/l2Tests/ (tests/output/l2Tests-poc/ when nReps = 1).


# ---------------------------------------------------------------------------
# Run settings: edit these
# ---------------------------------------------------------------------------

nReps <- 1          # replicates per setting: 1 = proof of concept, 50 for the paper
nCores <- parallel::detectCores() - 1

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

# sections <- as.integer(commandArgs(trailingOnly = TRUE))
# if (length(sections) == 0L) sections <- 0:3
sections <- 3

nReps <- as.integer(Sys.getenv("L2_TESTS_REPS", nReps))
nCores <- max(1L, as.integer(Sys.getenv("L2_TESTS_CORES", nCores)))

outputDir <- file.path("tests", "output", if (nReps == 1) "l2Tests-poc" else "l2Tests")
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

methodColours <- c("THMM (funHMM)" = "#000000", "THSMM, geometric" = "#0072B2",
                   "THSMM, shifted Poisson" = "#CC79A7", "THSMM, negative binomial" = "#D55E00",
                   "THSMM, nonparametric" = "#E69F00", "THSMM, true family" = "#009E73")   # Okabe-Ito
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
        for (m in methods) lines(seq_len(maxPlot), quantiles[[m]][2, ], type = "b", pch = 19, cex = 0.5, col = methodColours[[m]])
        if (i == 1) legend("topright", c("truth", methods), col = c("grey80", methodColours[methods]),
                           lwd = c(7, rep(1.5, length(methods))), bty = "n", cex = 0.75)
      }
      mtext(sprintf("%s: %s dwell times, %s separation (a = %g)%s", title, tolower(truthName), separation,
                    separations[[separation]], if (nReps > 1) sprintf(", median and 10-90%% band over %d replicates", nReps) else ""),
            outer = TRUE, font = 2)
    }
  }
}

cat(sprintf("l2Tests: sections %s, %d replicates, %d cores\n", paste(sections, collapse = ", "), nReps, nCores))


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
}


# ===========================================================================
# Section 1: recovering a geometric sojourn distribution
# ===========================================================================

if (1 %in% sections) {
  cat("\nSection 1: geometric sojourn recovery\n")
  result1 <- runTasks("Geometric", 1, c("THMM (funHMM)", "THSMM, geometric", "THSMM, negative binomial"))
  saveRDS(result1, file.path(outputDir, "section1.rds"))
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
  result2 <- runTasks(names(dwellTruths), 2, c("THMM (funHMM)", fittedFamilies))
  saveRDS(result2, file.path(outputDir, "section2.rds"))
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
  result3 <- runTasks(names(dwellTruths), 3, c("THMM (funHMM)", "THSMM, geometric", "THSMM, true family"))
  saveRDS(result3, file.path(outputDir, "section3.rds"))
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

cat(sprintf("\nResults written to %s\n", outputDir))
