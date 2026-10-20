# THMM vs THSMM on the EMG (EMG31) power spectra of two patients (CF050, CF076)
#
# Run from anywhere:   Rscript osaModels/emgModels.R
# Environment variable OSA_CORES sets the parallel workers (default detectCores() - 1).
#
# The same procedure as eegModelsRefined.R, applied to the chin EMG instead of
# the C4-M1 EEG: the same preprocessing (stage sequences for 2 / 3 / 5 states,
# log-power curve per 30 s epoch cut to Freq 0-49, ksmooth(bandwidth = 5) smoothing),
# the same settings and the same models (108 fits):
#   patient   CF050 and CF076 (data/Wpsd-raw-EMG31-*.csv; the stage scoring is
#             the same as in the EEG files)
#   nStates   2 (awake / asleep), 3 (awake / REM / NREM), 5 (scored stages)
#   pipeline  "smoothed W22" and "smoothed W21" on Freq 0-49, funHMM's data-driven
#             sigma converted to the package's emission scale (the pipelines
#             chosen from the EEG grid in eegModels.R)
#   models    THMM (funHMM), and THSMMs with Weibull, neg. binomial, Poisson and
#             npg(25) dwell, each from two unsupervised starts:
#               "k-means"  k-means mean curves, and transitions / mean dwells from
#                          the k-means labels (the THMM starts here too);
#               "THMM"     the fitted THMM: its mean curves, embedded transitions,
#                          mean dwells 1 / (1 - A_jj) and initial distribution, so
#                          the THSMM only has to improve on the THMM's solution.
#
# Scores of each Viterbi path against the scored stages: accuracy and balanced
# accuracy (mean per-stage recall) after the best matching of fitted states to
# stages, and the ARI. Log-likelihoods and BICs are comparable within a
# (patient, nStates, pipeline) setting only.
#
# Output in osaModels/outputEmg/: results.csv, accuracy / BIC / winner tables
# (.csv and LaTeX .tex), confusion.txt, emgModels.pdf, and a cache of finished
# fits (cache/; delete it after changing the EM settings).

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

if (!requireNamespace("funHMM", quietly = TRUE)) {
  stop("funHMM is not installed: install.packages(\"funHMM\")")
}
library(parallel)

nCores <- as.integer(Sys.getenv("OSA_CORES", max(1L, detectCores() - 1L)))
outputDir <- file.path("osaModels", "outputEmg")
cacheDir <- file.path(outputDir, "cache")
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Data: the preprocessing of eegModels.R, for each patient
# ---------------------------------------------------------------------------

nFrequencies <- 50   # Freq 0-49; cut before smoothing, so higher frequencies never enter

loadPatient <- function(file) {
  eeg <- read.csv(file)
  # True sleep stages (1 awake, 2 REM, 3-5 NREM 1-3) and the PSD curve of each epoch
  ss <- eeg$Stage[seq(1, length(eeg$Stage), 257)]
  dat <- matrix(eeg$Power, ncol = 257, byrow = TRUE)
  # Keep only the first nFrequencies frequencies (Freq 0-49) of each epoch
  dat <- dat[, seq_len(nFrequencies)]
  # 3 states: Awake, REM, Non-REM; 2 states: Awake, Asleep
  ss3 <- ss
  ss3[which(ss >= 3)] <- 3
  ss2 <- ss
  ss2[which(ss >= 2)] <- 2
  # Smoothed curves
  sm <- dat
  for (i in 1:nrow(sm)) {
    sm[i, ] <- ksmooth(seq_len(nFrequencies), sm[i, ], bandwidth = 5, x.points = seq_len(nFrequencies))$y
  }
  list(truths = list(`2` = ss2, `3` = ss3, `5` = ss), curves = list(raw = dat, smoothed = sm))
}

patients <- list(CF050 = loadPatient("./data/Wpsd-raw-EMG31-CF050.csv"),
                 CF076 = loadPatient("./data/Wpsd-raw-EMG31-CF076.csv"))


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

stageLabels <- list(`2` = c("Awake", "Asleep"), `3` = c("Awake", "REM", "NREM"),
                    `5` = c("Awake", "REM", "NREM1", "NREM2", "NREM3"))
nStatesValues <- c(2L, 3L, 5L)
pipelines <- list(`smoothed W22` = list(data = "smoothed", norm = "W22"),
                  `smoothed W21` = list(data = "smoothed", norm = "W21"))

parametricDwellSize <- 150        # exact dwell support; scored sojourns reach 169 epochs
seed <- 137
emTol <- 1e-8
emMaxIter <- 1000

thmmMethod <- "THMM (funHMM)"
thsmmFamilies <- list(
  Weibull         = list(type = "weibull", size = parametricDwellSize),
  `neg. binomial` = list(type = "nb", size = parametricDwellSize),
  Poisson         = list(type = "pois", size = parametricDwellSize),
  `npg(25)`       = list(type = "npg", size = 25)
)
startTypes <- c("k-means", "THMM")
thsmmMethods <- as.vector(outer(names(thsmmFamilies), startTypes,
                                function(f, s) sprintf("THSMM, %s (%s start)", f, s)))
methodOrder <- c(thmmMethod, thsmmMethods)
settingColumns <- c("patient", "nStates", "pipeline")


# ---------------------------------------------------------------------------
# Starting values and fits
# ---------------------------------------------------------------------------

# funHMM writes the functional as -||T(O - h)||^2 / (2 sigma^2) with T = 1 / sqrt(len)
# (L2), sqrt(len) * first differences (W21) or len^1.5 * second differences (W22);
# the package as -scale / 2 times the plain sum of squares.
normFactor <- function(norm, len) switch(norm, L2 = 1 / len, W21 = len, W22 = len^3)

emissionSettings <- function(x, norm) {
  sigma <- funHMM:::.sigma_auto(array(x, c(nrow(x), ncol(x), 1L)), "nonpar", norm, 0.5)
  list(sigma = sigma, control = list(norm = norm, scale = normFactor(norm, ncol(x)) / sigma^2))
}

kmeansStart <- function(x, nStates, control) {
  family <- getEmissionFamily("nonpar")
  set.seed(seed)
  init <- family$initialize(family$prepare(x, control), nStates, control, "kmeans")
  delta <- tabulate(init$labels[1], nbins = nStates) + 0.1
  list(means = init$params, tpm = initializeEmbeddedTpmFromLabels(init$labels, nStates),
       meanDwells = initializeMeanDwellsFromLabels(init$labels, nStates), delta = delta / sum(delta))
}

# The fitted THMM as a THSMM start: A = diag(stay) + (1 - stay) * tpm.
thmmStart <- function(thmm, nStates) {
  stay <- pmin(diag(thmm$trans), 1 - 1e-6)
  tpm <- thmm$trans
  diag(tpm) <- 0
  tpm <- if (nStates == 2L) matrix(c(0, 1, 1, 0), 2) else tpm / rowSums(tpm)
  list(means = thmm$means, tpm = tpm, meanDwells = 1 / (1 - stay), delta = thmm$delta)
}

fitFunHmm <- function(x, nStates, norm, sigma, start) {
  stay <- 1 - pmin(pmax(1 / start$meanDwells, 1e-3), 1 - 1e-3)   # geometric dwell with the start's means
  trans <- diag(stay, nStates) + (1 - stay) * start$tpm
  fit <- suppressWarnings(funHMM::thmm(x, nStates, type = "nonpar", norm = norm, sigma = sigma,
                                       par = start$means, trans = trans, init.prob = start$delta,
                                       tol = emTol, max.iter = emMaxIter))
  list(states = fit$states, logLik = fit$loglik, converged = fit$converged, iterations = fit$iter,
       nParameters = length(start$means) + nStates * (nStates - 1) + (nStates - 1),
       means = fit$par, trans = fit$trans, delta = fit$init.prob)
}

fitThsmm <- function(x, nStates, control, start, family) {
  sojournParameters <- lapply(seq_len(nStates), function(i) {
    initializeSojournOne(family$type, dwellSize = family$size, meanDwell = start$meanDwells[i])
  })
  model <- buildHsmm(start$tpm, sojournParameters, rep(family$type, nStates), start$means, "nonpar",
                     delta = start$delta, stationary = FALSE, emissionControl = control)
  fit <- suppressWarnings(fitHsmmEM(model, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
  fitted <- fit$fittedModel
  dwellCount <- sum(mapply(getSojournParamCount, fitted$sojournParameters, fitted$sojournType))
  list(states = hsmmViterbiEM(fitted, x)$states, logLik = fit$logLik, converged = fit$converged,
       iterations = fit$iterations,
       nParameters = length(fitted$emissionParameters) + dwellCount + nStates * max(nStates - 2, 0) + (nStates - 1),
       model = fitted[c("tpm", "sojournType", "sojournParameters", "emissionParameters", "delta")])
}


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

allPermutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  smaller <- allPermutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))
  }))
}
permutations <- lapply(setNames(nStatesValues, nStatesValues), allPermutations)

# Fitted states relabelled to the stages they match best (maximum agreement).
alignToTruth <- function(states, truth, nStates) {
  tab <- table(factor(states, seq_len(nStates)), factor(truth, seq_len(nStates)))
  perms <- permutations[[as.character(nStates)]]
  best <- perms[which.max(apply(perms, 1, function(p) sum(tab[cbind(p, seq_len(nStates))]))), ]
  match(states, best)
}

scorePath <- function(states, truth, nStates) {
  aligned <- alignToTruth(states, truth, nStates)
  recall <- vapply(seq_len(nStates), function(k) mean(aligned[truth == k] == k), numeric(1))
  c(accuracy = mean(aligned == truth), balancedAccuracy = mean(recall),
    ARI = mclust::adjustedRandIndex(truth, states))
}


# ---------------------------------------------------------------------------
# Fits: one job per setting (the THMM-start THSMMs need that setting's THMM)
# ---------------------------------------------------------------------------

settings <- expand.grid(pipeline = names(pipelines), nStates = nStatesValues, patient = names(patients),
                        stringsAsFactors = FALSE)[, settingColumns]
settings <- settings[order(-settings$nStates), ]   # longest first
rownames(settings) <- NULL

cachedFit <- function(setting, method, fitFunction) {
  key <- paste(c(setting, method, nFrequencies, emTol, emMaxIter, parametricDwellSize), collapse = "_")
  file <- file.path(cacheDir, paste0(gsub("[^A-Za-z0-9._-]", "", gsub("[ ,()]+", "-", key)), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  elapsed <- system.time(fit <- tryCatch(fitFunction(), error = function(e) list(error = conditionMessage(e))))
  fit$seconds <- elapsed[["elapsed"]]
  saveRDS(fit, file)
  fit
}

runSetting <- function(setting) {
  pipeline <- pipelines[[setting$pipeline]]
  nStates <- setting$nStates
  x <- patients[[setting$patient]]$curves[[pipeline$data]]
  emission <- emissionSettings(x, pipeline$norm)
  starts <- list(`k-means` = kmeansStart(x, nStates, emission$control))
  fits <- list()
  fits[[thmmMethod]] <- cachedFit(setting, thmmMethod, function() {
    fitFunHmm(x, nStates, pipeline$norm, emission$sigma, starts[["k-means"]])
  })
  if (is.null(fits[[thmmMethod]]$error)) starts$THMM <- thmmStart(fits[[thmmMethod]], nStates)
  for (startType in startTypes) {
    for (family in names(thsmmFamilies)) {
      method <- sprintf("THSMM, %s (%s start)", family, startType)
      fits[[method]] <- if (is.null(starts[[startType]])) list(error = "no THMM start") else
        cachedFit(setting, method, function() {
          fitThsmm(x, nStates, emission$control, starts[[startType]], thsmmFamilies[[family]])
        })
    }
  }
  fits
}

cat(sprintf("emgModels: %d settings, %d fits on %d cores\n", nrow(settings),
            nrow(settings) * length(methodOrder), nCores))
started <- Sys.time()
fits <- mclapply(seq_len(nrow(settings)), function(s) {
  tryCatch(runSetting(settings[s, ]), error = function(e) list(error = conditionMessage(e)))
}, mc.cores = nCores, mc.preschedule = FALSE)
cat(sprintf("Fitted in %.1f minutes\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

isFailed <- function(fit) !is.list(fit) || !is.null(fit$error) || is.null(fit$states) || !is.finite(fit$logLik)
nEpochs <- vapply(patients, function(p) length(p$truths[[1]]), integer(1))

results <- do.call(rbind, lapply(seq_len(nrow(settings)), function(s) {
  setting <- settings[s, ]
  nStates <- setting$nStates
  truth <- patients[[setting$patient]]$truths[[as.character(nStates)]]
  do.call(rbind, lapply(methodOrder, function(method) {
    fit <- fits[[s]][[method]] %||% list(error = fits[[s]]$error %||% "not fitted")
    failed <- isFailed(fit)
    score <- if (failed) c(accuracy = NA, balancedAccuracy = NA, ARI = NA) else scorePath(fit$states, truth, nStates)
    data.frame(setting, method = method,
               family = if (method == thmmMethod) "THMM" else sub("^THSMM, (.*) \\(.*$", "\\1", method),
               start = if (method == thmmMethod) "k-means" else sub("^.*\\((.*) start\\)$", "\\1", method),
               t(score),
               logLik = if (failed) NA else fit$logLik,
               nParameters = if (failed) NA else fit$nParameters,
               BIC = if (failed) NA else -2 * fit$logLik + fit$nParameters * log(nEpochs[[setting$patient]]),
               converged = if (failed) NA else fit$converged,
               iterations = if (failed) NA else fit$iterations,
               seconds = round(if (is.list(fit)) fit$seconds %||% NA else NA, 1),
               error = if (!is.list(fit)) as.character(fit) else fit$error %||% "",
               row.names = NULL)
  }))
}))
results <- results[order(results$patient, results$nStates, results$pipeline, match(results$method, methodOrder)), ]
results$setting <- paste(results$patient, results$nStates, results$pipeline, sep = " | ")

# Within each setting: rank by accuracy (1 = best, ties share the rank) and gain over the THMM.
results$rank <- NA_real_
results$gainOverThmm <- NA_real_
for (rows in split(seq_len(nrow(results)), results$setting)) {
  results$rank[rows] <- rank(-results$accuracy[rows], na.last = "keep", ties.method = "min")
  thmm <- results$accuracy[rows][results$method[rows] == thmmMethod]
  if (length(thmm) == 1L) results$gainOverThmm[rows] <- results$accuracy[rows] - thmm
}
rownames(results) <- NULL

failures <- results[results$error != "" | is.na(results$accuracy), ]
if (nrow(failures)) {
  cat(sprintf("%d of %d fits failed:\n", nrow(failures), nrow(results)))
  print(failures[, c("setting", "method", "error")], row.names = FALSE)
}
unconverged <- sum(!results$converged, na.rm = TRUE)
if (unconverged) cat(sprintf("%d fits reached emMaxIter = %d without converging\n", unconverged, emMaxIter))

saveRDS(list(settings = settings, fits = fits, results = results), file.path(outputDir, "fits.rds"))
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

# One row per setting, one column per model.
wideTable <- function(column) {
  wide <- reshape(results[, c(settingColumns, "method", column)], idvar = settingColumns,
                  timevar = "method", direction = "wide")
  names(wide) <- sub(paste0("^", column, "\\."), "", names(wide))
  wide <- wide[, c(settingColumns, intersect(methodOrder, names(wide)))]
  names(wide) <- sub("^THSMM, ", "", names(wide))
  rownames(wide) <- NULL
  wide
}

# Per setting: the THMM, the most accurate THSMM (ties broken by BIC), and the THSMM with the smallest BIC.
settingBest <- do.call(rbind, lapply(split(results, results$setting), function(r) {
  thsmm <- r[r$method != thmmMethod & !is.na(r$accuracy), ]
  thmm <- r$accuracy[r$method == thmmMethod]
  if (!nrow(thsmm)) return(NULL)
  best <- thsmm[order(-thsmm$accuracy, thsmm$BIC), ][1, ]
  bic <- thsmm[which.min(thsmm$BIC), ]
  data.frame(r[1, settingColumns], thmmAccuracy = thmm, bestThsmm = best$method, bestAccuracy = best$accuracy,
             gain = best$accuracy - thmm, bicThsmm = bic$method, bicAccuracy = bic$accuracy,
             bicGain = bic$accuracy - thmm, row.names = NULL)
}))
settingBest <- settingBest[order(settingBest$patient, settingBest$nStates, settingBest$pipeline), ]
rownames(settingBest) <- NULL

# Per patient and over both: mean accuracy, mean rank, and how often / by how much each THSMM beats the THMM.
summarizeMethods <- function(rows, label) {
  out <- do.call(rbind, lapply(split(rows, rows$method), function(g) {
    data.frame(patient = label, method = g$method[1], settings = nrow(g), failed = sum(is.na(g$accuracy)),
               accuracy = mean(g$accuracy, na.rm = TRUE), balancedAccuracy = mean(g$balancedAccuracy, na.rm = TRUE),
               ARI = mean(g$ARI, na.rm = TRUE), meanRank = mean(g$rank, na.rm = TRUE),
               gainOverThmm = mean(g$gainOverThmm, na.rm = TRUE), beatsThmm = sum(g$gainOverThmm > 0, na.rm = TRUE),
               losesToThmm = sum(g$gainOverThmm < 0, na.rm = TRUE), row.names = NULL)
  }))
  out[order(out$meanRank), ]
}
methodSummary <- rbind(do.call(rbind, lapply(names(patients), function(p) summarizeMethods(results[results$patient == p, ], p))),
                       summarizeMethods(results, "both"))
rownames(methodSummary) <- NULL

tables <- list(
  accuracy = list(table = wideTable("accuracy"), caption = paste(
    "Sleep-stage EMG (EMG31): Viterbi accuracy against the scored stages after matching states to stages.")),
  balancedAccuracy = list(table = wideTable("balancedAccuracy"), caption = paste(
    "Sleep-stage EMG: balanced accuracy (mean per-stage recall) of the Viterbi path.")),
  bic = list(table = wideTable("BIC"), caption = paste(
    "Sleep-stage EMG: BIC under the nonparametric emission functional; comparable within a row only.")),
  settingBest = list(table = settingBest, caption = paste(
    "Per setting: funHMM's THMM, the most accurate THSMM, and the THSMM with the smallest BIC.")),
  methodSummary = list(table = methodSummary, caption = paste(
    "Each model averaged over the settings of each patient: accuracy, mean accuracy rank within a setting,",
    "mean gain over the THMM, and the number of settings in which it beats / loses to the THMM."))
)
for (name in names(tables)) {
  write.csv(tables[[name]]$table, file.path(outputDir, paste0(name, ".csv")), row.names = FALSE)
  writeLatex(tables[[name]]$table, file.path(outputDir, paste0(name, ".tex")), tables[[name]]$caption)
}

printTable <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(format(table, digits = 3), row.names = FALSE)
}
printTable("Viterbi accuracy", tables$accuracy$table)
printTable("THMM vs the best and the BIC-chosen THSMM", settingBest)
printTable("Models averaged over settings", methodSummary)
for (p in c(names(patients), "both")) {
  winner <- methodSummary[methodSummary$patient == p & methodSummary$method != thmmMethod, ][1, ]
  cat(sprintf("Winning THSMM, %s: %s (mean rank %.2f, beats the THMM in %d of %d settings, mean gain %+.3f)\n",
              p, winner$method, winner$meanRank, winner$beatsThmm, winner$settings, winner$gainOverThmm))
}

sink(file.path(outputDir, "confusion.txt"))
for (s in seq_len(nrow(settingBest))) {
  setting <- settingBest[s, ]
  nStates <- setting$nStates
  truth <- patients[[setting$patient]]$truths[[as.character(nStates)]]
  fitList <- fits[[which(settings$patient == setting$patient & settings$nStates == nStates &
                           settings$pipeline == setting$pipeline)]]
  for (method in c(thmmMethod, setting$bestThsmm)) {
    if (isFailed(fitList[[method]])) next
    cat(sprintf("\n%s | %d states | %s | %s\n", setting$patient, nStates, setting$pipeline, method))
    print(table(stage = factor(truth, seq_len(nStates), stageLabels[[as.character(nStates)]]),
                viterbi = factor(alignToTruth(fitList[[method]]$states, truth, nStates), seq_len(nStates),
                                 stageLabels[[as.character(nStates)]])))
  }
}
sink()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

familyColours <- c(THMM = "#000000", Weibull = "#CC79A7", `neg. binomial` = "#D55E00", Poisson = "#E69F00",
                   `npg(25)` = "#009E73")

pdf(file.path(outputDir, "emgModels.pdf"), width = 11, height = 7)

# Accuracy of every model in every setting, one panel per patient.
par(mfrow = c(1, length(patients)), mar = c(9, 4, 3, 1), oma = c(0, 0, 2, 0))
for (p in names(patients)) {
  table <- tables$accuracy$table[tables$accuracy$table$patient == p, ]
  values <- t(as.matrix(table[, sub("^THSMM, ", "", methodOrder)]))
  methodFamily <- results$family[match(methodOrder, results$method)]
  colours <- familyColours[methodFamily]
  density <- ifelse(grepl("THMM start", methodOrder), 30, -1)
  barplot(values, beside = TRUE, col = colours, density = density, names.arg = paste(table$nStates, "st.,", table$pipeline),
          las = 2, ylim = c(0, 1), ylab = "Viterbi accuracy", main = p)
  if (p == names(patients)[1]) {
    legend("topright", c(names(familyColours), "hatched: THMM start"), fill = c(familyColours, "grey"),
           density = c(rep(-1, length(familyColours)), 30), bty = "n", cex = 0.7)
  }
}
mtext("Viterbi accuracy by setting and model", outer = TRUE, font = 2)

# Hypnograms: scored stages, the THMM and the best THSMM, for every patient and nStates
# in the pipeline where that THSMM is most accurate.
for (p in names(patients)) {
  for (nStates in nStatesValues) {
    rows <- settingBest[settingBest$patient == p & settingBest$nStates == nStates, ]
    if (!nrow(rows)) next
    setting <- rows[which.max(rows$bestAccuracy), ]
    truth <- patients[[p]]$truths[[as.character(nStates)]]
    fitList <- fits[[which(settings$patient == p & settings$nStates == nStates & settings$pipeline == setting$pipeline)]]
    labels <- stageLabels[[as.character(nStates)]]
    par(mfrow = c(3, 1), mar = c(1, 6, 1.5, 1), oma = c(3, 0, 2, 0))
    plot(truth, type = "s", yaxt = "n", xlab = "", ylab = "", main = "Scored stages", ylim = c(0.5, nStates + 0.5))
    axis(2, at = seq_len(nStates), labels = labels, las = 1, cex.axis = 0.7)
    for (method in c(thmmMethod, setting$bestThsmm)) {
      fit <- fitList[[method]]
      if (isFailed(fit)) { plot.new(); next }
      accuracy <- results$accuracy[results$setting == paste(p, nStates, setting$pipeline, sep = " | ") &
                                     results$method == method]
      colour <- familyColours[[results$family[match(method, results$method)]]]
      plot(alignToTruth(fit$states, truth, nStates), type = "s", col = colour, yaxt = "n", xlab = "", ylab = "",
           ylim = c(0.5, nStates + 0.5), main = sprintf("%s (accuracy %.3f)", method, accuracy))
      axis(2, at = seq_len(nStates), labels = labels, las = 1, cex.axis = 0.7)
    }
    mtext(sprintf("%s, %d states, %s", p, nStates, setting$pipeline), outer = TRUE, font = 2)
    mtext("epoch", side = 1, outer = TRUE, line = 1.5)
  }
}

# Fitted dwell pmfs of the best 5-state THSMM against the scored sojourn lengths.
for (p in names(patients)) {
  rows <- settingBest[settingBest$patient == p & settingBest$nStates == 5L, ]
  if (!nrow(rows)) next
  setting <- rows[which.max(rows$bestAccuracy), ]
  fit <- fits[[which(settings$patient == p & settings$nStates == 5L & settings$pipeline == setting$pipeline)]][[setting$bestThsmm]]
  if (isFailed(fit)) next
  truth <- patients[[p]]$truths[["5"]]
  runs <- rle(truth)
  aligned <- alignToTruth(fit$states, truth, 5L)
  model <- buildHsmm(fit$model$tpm, fit$model$sojournParameters, fit$model$sojournType, fit$model$emissionParameters,
                     "nonpar", delta = fit$model$delta, stationary = FALSE,
                     emissionControl = list(norm = pipelines[[setting$pipeline]]$norm))
  colour <- familyColours[[results$family[match(setting$bestThsmm, results$method)]]]
  par(mfrow = c(1, 5), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (s in 1:5) {
    observed <- runs$lengths[runs$values == s]
    breaks <- seq(0, max(observed, 10) + 10, by = 5)
    hist(observed, breaks = breaks, freq = FALSE, col = "grey85", border = "white", xlab = "sojourn (epochs)",
         main = sprintf("%s (n = %d)", stageLabels[["5"]][s], length(observed)))
    fittedState <- unique(fit$states[aligned == s])
    if (length(fittedState)) {
      pmf <- getModelDwellPmf(model, fittedState[1], max(breaks))
      lines(seq_along(pmf), pmf, col = colour, lwd = 2)
    }
  }
  mtext(sprintf("%s, 5 states, %s: scored sojourns and the dwell pmf of %s", p, setting$pipeline, setting$bestThsmm),
        outer = TRUE, font = 2)
}

invisible(dev.off())
cat(sprintf("\nResults written to %s\n", outputDir))
