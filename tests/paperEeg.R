# Sleep-stage EEG application: THMM vs THSMM on the C4-M1 power spectra of
# subject CF050 (data/Wpsd-raw-EEGC4M1-CF050.csv).
#
# Run from the package root:   Rscript tests/paperEeg.R
# Environment variable PAPER_EEG_CORES sets the parallel workers (default
# detectCores() - 1).
#
# The THMM fits follow archivedCode/thmmEEG.r exactly:
#   * one curve per 30 s epoch: the 257-point log-power spectrum, of which
#     only the first 50 frequencies are used;
#   * raw curves, and curves smoothed by ksmooth(bandwidth = 5);
#   * nonparametric state means with the L2, W21 and W22 norms;
#   * 2 states (truth: awake vs asleep) with prob.trans = 0.8 I + 0.1, tol 1e-12
#     and thmm.BaumWelch()'s default maxIter = 100; and 5 states (truth: the
#     five scored stages) with prob.trans = 0.5 I + 0.1, tol 1e-13, maxIter 1000;
#   * starting means from autoParam() (random epochs) and a Dirichlet initial
#     distribution, both drawn after set.seed(137). (thmmEEG.r sets the seed
#     only for its 5-state fits; here every fit is seeded, so all runs are
#     reproducible.) The starts are drawn exactly as thmm.BaumWelch() draws
#     them, so passing them in reproduces thmmEEG.r's own call.
#
# Every model below starts from those same mean curves, initial distribution
# and transition structure (embedded transitions omega and mean dwell
# 1 / (1 - A_ii) taken from prob.trans), and is scored by the ARI and the
# cross-entropy of its Viterbi path against the scored stages:
#   THMM (thmm1.R)          thmmEEG.r's fit, thmm1.R unmodified
#   THMM (package)          geometric dwell, fitted by EM to convergence. Same
#                           functional as thmm1.R, so identical iterates
#                           (tests/test-geo-vs-thmm1.R); only the stopping rule differs
#   THSMM, npg(10)          nonparametric dwell head + geometric tail
#   THSMM, Poisson          1 + Poisson dwell
#   THSMM, neg. binomial    1 + negative binomial dwell
#   THSMM, Weibull          discretized Weibull dwell
#   THSMM, shifted Poisson  shift profiled over shiftCandidates by
#                           fitHsmmShiftGrid (2-state models only: the 5-state
#                           grid has too many combinations)
# Log-likelihoods are those of the nonparametric emission functional, so they
# (and the BIC) compare models with the same norm and data only.
#
# Output in tests/output/paperEeg/: ARI / cross-entropy / BIC tables (.csv and
# LaTeX .tex), confusion tables (confusion.txt), all fits (fits.rds) and
# figures (paperEeg.pdf).

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

thmmPath <- file.path("..", "archivedCode", "thmm1.r")
if (!file.exists(thmmPath)) thmmPath <- file.path("tests", "thmm1.r")
if (!file.exists(thmmPath)) stop("thmm1.r not found in ../archivedCode or tests/")
thmm1 <- new.env()
suppressMessages(source(thmmPath, local = thmm1))   # also attaches mclust

library(parallel)
nCores <- as.integer(Sys.getenv("PAPER_EEG_CORES", max(1L, detectCores() - 1L)))

outputDir <- file.path("tests", "output", "paperEeg")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Data, as in thmmEEG.r
# ---------------------------------------------------------------------------

nFrequencies <- 257
usedFrequencies <- 1:50

eeg <- read.csv(file.path("data", "Wpsd-raw-EEGC4M1-CF050.csv"))
stages <- eeg$Stage[seq(1, length(eeg$Stage), nFrequencies)]   # 1 awake, 2 REM, 3-5 NREM 1-3
curves <- matrix(eeg$Power, ncol = nFrequencies, byrow = TRUE)
smoothed <- curves
for (i in seq_len(nrow(smoothed))) {
  smoothed[i, ] <- ksmooth(seq_len(nFrequencies), smoothed[i, ], bandwidth = 5)$y
}

datasets <- list(raw = curves[, usedFrequencies], smoothed = smoothed[, usedFrequencies])
truths <- list(`2` = pmin(stages, 2L), `5` = stages)   # awake / asleep, and the five stages
stageLabels <- list(`2` = c("Awake", "Asleep"), `5` = c("Awake", "REM", "NREM1", "NREM2", "NREM3"))

thmmSettings <- list(
  `2` = list(prob.trans = 0.8 * diag(2) + matrix(0.1, 2, 2), tol = 1e-12, maxIter = 100),
  `5` = list(prob.trans = 0.5 * diag(5) + matrix(0.1, 5, 5), tol = 1e-13, maxIter = 1000)
)

norms <- c("L2", "W21", "W22")
seed <- 137

npgSize <- 10
parametricDwellSize <- 150        # exact dwell support; observed sojourns reach 355 epochs
shiftCandidates <- 1:6
emTol <- 1e-10
emMaxIter <- 2000

thsmmFamilies <- c("THSMM, npg(10)" = "npg", "THSMM, Poisson" = "pois", "THSMM, neg. binomial" = "nb",
                   "THSMM, Weibull" = "weibull", "THSMM, shifted Poisson" = "shiftp")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Cross-entropy of articleSims.r, called as in thmmEEG.r: crossEnt(truth, estimate).
crossEnt <- function(est, tru) {
  tab <- table(est, tru)
  ftab <- t(t(tab) / colSums(tab))
  -sum(apply(log2(ftab), 2, max))
}

# thmm.BaumWelch()'s own random start (autoParam, then thmm.init's Dirichlet
# initial distribution), drawn explicitly so every model can share it.
drawStart <- function(x, nStates) {
  set.seed(seed)
  means <- thmm1$autoParam(x, nStates, "nonpar")
  delta <- as.vector(thmm1$rdirichlet(1, rep(6, nStates)))
  list(means = means, delta = delta)
}

fitThmm1 <- function(x, nStates, norm, start) {
  setting <- thmmSettings[[as.character(nStates)]]
  invisible(capture.output(fit <- thmm1$thmm.BaumWelch(
    data = x, stateCnt = nStates, state_mean = start$means, prob.trans = setting$prob.trans,
    prob.init = start$delta, tol = setting$tol, maxIter = setting$maxIter, verbose = 1,
    type = "nonpar", PLOT = FALSE, norm = norm
  )))
  # thmm.BaumWelch() stores alfa from before its last re-estimation, so the
  # log-likelihood is recomputed at the returned parameters.
  final <- thmm1$thmm.alfa(thmm1$thmm.bj.all(fit, log.scale = TRUE, type = "nonpar", norm = norm), TRUE, "nonpar")
  list(states = fit$bestSeq, logLik = thmm1$sumLogs.vec(final$alfa[nrow(final$alfa), ]),
       nParameters = length(fit$stPar) + nStates * (nStates - 1) + (nStates - 1),
       A = exp(fit$prTran), means = fit$stPar)
}

# Package starting model with the THMM start: omega and mean dwells from prob.trans.
startModel <- function(nStates, norm, start, sojournType) {
  A <- thmmSettings[[as.character(nStates)]]$prob.trans
  stay <- diag(A)
  omega <- A / (1 - stay)
  diag(omega) <- 0
  dwellSize <- if (sojournType == "npg") npgSize else parametricDwellSize
  sojournParameters <- lapply(seq_len(nStates), function(i) {
    if (sojournType == "geo") return(c(prob = 1 - stay[i]))
    initializeSojournOne(sojournType, dwellSize = dwellSize, meanDwell = 1 / (1 - stay[i]))
  })
  buildHsmm(omega, sojournParameters, rep(sojournType, nStates), start$means, "nonpar",
            delta = start$delta, stationary = FALSE, emissionControl = list(norm = norm))
}

fitPackage <- function(x, nStates, norm, start, sojournType) {
  model <- startModel(nStates, norm, start, sojournType)
  fit <- if (sojournType == "shiftp") {
    suppressWarnings(fitHsmmShiftGrid(model, x, rep(list(shiftCandidates), nStates), stationary = FALSE,
                                      fitFunction = fitHsmmEM, tol = emTol, maxIter = emMaxIter))
  } else {
    suppressWarnings(fitHsmmEM(model, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
  }
  fitted <- fit$fittedModel
  dwellCount <- sum(mapply(getSojournParamCount, fitted$sojournParameters, fitted$sojournType,
                           MoreArgs = list(fitShift = FALSE)))
  list(states = hsmmViterbiEM(fitted, x)$states, logLik = fit$logLik,
       nParameters = length(fitted$emissionParameters) + dwellCount + nStates * max(nStates - 2, 0) +
         (nStates - 1) + if (sojournType == "shiftp") nStates else 0,   # shifts are profiled
       model = fitted, converged = fit$converged, iterations = fit$iterations, bestShift = fit$bestShift)
}


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

methodOrder <- c("THMM (thmm1.R)", "THMM (package)", names(thsmmFamilies))
jobs <- expand.grid(method = methodOrder,
                    norm = norms, data = names(datasets), nStates = c(2L, 5L), stringsAsFactors = FALSE)
jobs <- jobs[!(jobs$method == "THSMM, shifted Poisson" & jobs$nStates == 5L), ]
jobs <- jobs[order(jobs$method == "THSMM, shifted Poisson", decreasing = TRUE), ]   # longest first

cat(sprintf("paperEeg: %d epochs, %d fits on %d cores\n", nrow(curves), nrow(jobs), nCores))
started <- Sys.time()

fits <- mclapply(seq_len(nrow(jobs)), function(j) {
  job <- jobs[j, ]
  x <- datasets[[job$data]]
  start <- drawStart(x, job$nStates)
  elapsed <- system.time(fit <- tryCatch({
    if (job$method == "THMM (thmm1.R)") {
      fitThmm1(x, job$nStates, job$norm, start)
    } else {
      sojournType <- if (job$method == "THMM (package)") "geo" else thsmmFamilies[[job$method]]
      fitPackage(x, job$nStates, job$norm, start, sojournType)
    }
  }, error = function(e) list(error = conditionMessage(e))))[["elapsed"]]
  fit$seconds <- elapsed
  fit
}, mc.cores = nCores, mc.preschedule = FALSE)

cat(sprintf("Fitted in %.1f minutes\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))

results <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  fit <- fits[[j]]
  truth <- truths[[as.character(jobs$nStates[j])]]
  failed <- !is.list(fit) || !is.null(fit$error) || is.null(fit$states)
  if (failed) cat(sprintf("  FAILED: %s, %s, %s, %d states: %s\n", jobs$method[j], jobs$data[j], jobs$norm[j],
                          jobs$nStates[j], if (is.list(fit)) fit$error else as.character(fit)))
  data.frame(
    jobs[j, c("nStates", "data", "norm", "method")],
    ARI = if (failed) NA else adjustedRandIndex(truth, fit$states),
    crossEnt = if (failed) NA else crossEnt(truth, fit$states),
    logLik = if (failed) NA else fit$logLik,
    nParameters = if (failed) NA else fit$nParameters,
    BIC = if (failed) NA else -2 * fit$logLik + fit$nParameters * log(nrow(curves)),
    seconds = round(fit$seconds %||% NA, 1),
    row.names = NULL
  )
}))
results <- results[order(results$nStates, results$data, results$norm, match(results$method, methodOrder)), ]
rownames(results) <- NULL

saveRDS(list(jobs = jobs, fits = fits, results = results), file.path(outputDir, "fits.rds"))
write.csv(results, file.path(outputDir, "results.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

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

wideTable <- function(column) {
  wide <- reshape(results[, c("nStates", "data", "norm", "method", column)],
                  idvar = c("nStates", "data", "norm"), timevar = "method", direction = "wide")
  names(wide) <- sub(paste0("^", column, "\\."), "", names(wide))
  wide <- wide[, c("nStates", "data", "norm", intersect(methodOrder, names(wide)))]
  rownames(wide) <- NULL
  wide
}

tables <- list(
  ari = list(table = wideTable("ARI"), caption = paste(
    "Sleep-stage EEG (subject CF050, C4-M1, first 50 frequencies): ARI of the Viterbi path against",
    "the scored stages (2 states: awake / asleep; 5 states: awake, REM, NREM1-3).")),
  crossEntropy = list(table = wideTable("crossEnt"), caption = paste(
    "Sleep-stage EEG: cross-entropy of the Viterbi path against the scored stages (lower is better).")),
  bic = list(table = wideTable("BIC"), caption = paste(
    "Sleep-stage EEG: BIC under the nonparametric emission functional; comparable within a row only."))
)

for (name in names(tables)) {
  write.csv(tables[[name]]$table, file.path(outputDir, paste0(name, ".csv")), row.names = FALSE)
  writeLatex(tables[[name]]$table, file.path(outputDir, paste0(name, ".tex")), tables[[name]]$caption)
}

cat("\nARI of the Viterbi path vs the scored stages\n")
print(format(tables$ari$table, digits = 3), row.names = FALSE)
cat("\nCross-entropy (lower is better)\n")
print(format(tables$crossEntropy$table, digits = 3), row.names = FALSE)
cat("\nBIC (within a row)\n")
print(format(tables$bic$table, digits = 6), row.names = FALSE)

# Confusion tables, as printed in thmmEEG.r: table(truth, Viterbi path).
sink(file.path(outputDir, "confusion.txt"))
for (j in seq_len(nrow(jobs))) {
  fit <- fits[[j]]
  if (!is.list(fit) || is.null(fit$states)) next
  truth <- truths[[as.character(jobs$nStates[j])]]
  cat(sprintf("\n%s | %s curves | %s | %d states%s\n", jobs$method[j], jobs$data[j], jobs$norm[j], jobs$nStates[j],
              if (!is.null(fit$bestShift)) paste(" | shifts", paste(fit$bestShift, collapse = ", ")) else ""))
  print(table(stage = factor(truth, seq_len(jobs$nStates[j]), stageLabels[[as.character(jobs$nStates[j])]]),
              viterbi = fit$states))
}
sink()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

methodColours <- c("THMM (thmm1.R)" = "#4d4d4d", "THMM (package)" = "#2166ac", "THSMM, npg(10)" = "#b2182b",
                   "THSMM, Poisson" = "#e08214", "THSMM, neg. binomial" = "#1a9850",
                   "THSMM, Weibull" = "#762a83", "THSMM, shifted Poisson" = "#35978f")

# Relabel a decoded path to the stages it overlaps most (for plotting only).
alignToTruth <- function(states, truth, nStates) {
  tab <- table(factor(states, seq_len(nStates)), factor(truth, seq_len(nStates)))
  permutations <- function(n) {
    if (n == 1L) return(matrix(1L, 1, 1))
    smaller <- permutations(n - 1L)
    do.call(rbind, lapply(seq_len(n), function(first) cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))))
  }
  all <- permutations(nStates)
  best <- all[which.max(apply(all, 1, function(p) sum(tab[cbind(p, seq_len(nStates))]))), ]
  match(states, best)
}

pdf(file.path(outputDir, "paperEeg.pdf"), width = 11, height = 7)

# ARI by setting and method.
par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
for (nStates in c(2L, 5L)) {
  table <- tables$ari$table[tables$ari$table$nStates == nStates, ]
  values <- t(as.matrix(table[, intersect(methodOrder, names(table))]))
  barplot(values, beside = TRUE, col = methodColours[rownames(values)], names.arg = paste(table$data, table$norm),
          las = 2, ylim = c(0, 1), ylab = "ARI", main = sprintf("%d states: Viterbi ARI", nStates))
  if (nStates == 2L) legend("topleft", rownames(values), fill = methodColours[rownames(values)], bty = "n", cex = 0.7)
}

# Decoded stage paths, as the last plot of thmmEEG.r, for the smoothed W22 fits.
for (nStates in c(2L, 5L)) {
  truth <- truths[[as.character(nStates)]]
  shown <- which(jobs$nStates == nStates & jobs$data == "smoothed" & jobs$norm == "W22")
  shown <- shown[order(match(jobs$method[shown], methodOrder))]
  par(mfrow = c(length(shown), 1), mar = c(1, 6, 1.5, 1), oma = c(3, 0, 2, 0))
  for (j in shown) {
    fit <- fits[[j]]
    plot(truth, type = "l", lwd = 2, yaxt = "n", xlab = "", ylab = "", ylim = c(0.5, nStates + 0.5),
         main = sprintf("%s (ARI %.3f)", jobs$method[j], results$ARI[results$method == jobs$method[j] &
                        results$nStates == nStates & results$data == "smoothed" & results$norm == "W22"]),
         cex.main = 0.9)
    if (is.list(fit) && !is.null(fit$states)) {
      lines(alignToTruth(fit$states, truth, nStates), col = methodColours[[jobs$method[j]]], lty = 2, lwd = 1.5)
    }
    axis(2, at = seq_len(nStates), labels = stageLabels[[as.character(nStates)]], las = 1, cex.axis = 0.7)
  }
  mtext(sprintf("%d-state models, smoothed curves, W22 norm: scored stages (black) and Viterbi path", nStates),
        outer = TRUE, font = 2)
  mtext("epoch", side = 1, outer = TRUE, line = 1.5)
}

# Fitted dwell distributions vs the observed sojourn lengths, 2-state smoothed W22.
truth <- truths[["2"]]
runs <- rle(truth)
shown <- which(jobs$nStates == 2L & jobs$data == "smoothed" & jobs$norm == "W22" & jobs$method != "THMM (thmm1.R)")
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (stage in 1:2) {
  observed <- runs$lengths[runs$values == stage]
  breaks <- seq(0, max(observed, 10) + 10, by = 10)
  hist(observed, breaks = breaks, freq = FALSE, col = "grey85", border = "white", xlab = "sojourn length (epochs)",
       main = sprintf("%s: observed sojourns (n = %d)", stageLabels[["2"]][stage], length(observed)))
  for (j in shown) {
    fit <- fits[[j]]
    if (!is.list(fit) || is.null(fit$model)) next
    # Fitted state that overlaps this stage most.
    tab <- table(factor(fit$states, 1:2), factor(truth, 1:2))
    fittedState <- if (tab[1, stage] + tab[2, 3 - stage] >= tab[2, stage] + tab[1, 3 - stage]) 1L else 2L
    pmf <- getModelDwellPmf(fit$model, fittedState, max(breaks))
    lines(seq_along(pmf), pmf, col = methodColours[[jobs$method[j]]], lwd = 1.5)
  }
  if (stage == 1) legend("topright", jobs$method[shown], col = methodColours[jobs$method[shown]], lwd = 1.5,
                         bty = "n", cex = 0.7)
}
mtext("2-state models, smoothed curves, W22 norm: fitted dwell pmfs", outer = TRUE, font = 2)

invisible(dev.off())

cat(sprintf("\nResults written to %s\n", outputDir))
