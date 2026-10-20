# THSMM osa model fitting and comparison 

# Load Source 
setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

if (!requireNamespace("funHMM", quietly = TRUE)) {
  stop("funHMM is not installed: install.packages(\"funHMM\")")
}

# Load in Data
eegc50 = read.csv("./data/Wpsd-raw-EEGC4M1-CF050.csv")

# Extract True Sleep State Sequences
eeg50.ss  = eegc50$Stage[seq(1,length(eegc50$Stage),257)]
# Extract PSD curvess for each Epoch
eeg50.dat = matrix(eegc50$Power,ncol = 257,byrow = T)
# Keep only the first 50 frequencies (Freq 0-49) of each epoch
eeg50.dat = eeg50.dat[,1:50]
# Compute the means for each sleep state
eeg50.mu  = rbind(
  apply(eeg50.dat[eeg50.ss==1,],2,mean),
  apply(eeg50.dat[eeg50.ss==2,],2,mean),
  apply(eeg50.dat[eeg50.ss==3,],2,mean),
  apply(eeg50.dat[eeg50.ss==4,],2,mean),
  apply(eeg50.dat[eeg50.ss==5,],2,mean)
)
# Create a 3 state sequence: Awake, REM, Non-REM
eeg50.ss3 <- eeg50.ss;
eeg50.ss3[which(eeg50.ss>=3)] <-3
# Create a 2 state sequence: Awake, Asleep
eeg50.ss2 <- eeg50.ss;
eeg50.ss2[which(eeg50.ss>=2)] <-2
# Smooth data
eeg50.sm = eeg50.dat
for( i in 1:nrow(eeg50.sm) ){
  eeg50.sm[i,] = ksmooth(1:ncol(eeg50.sm),eeg50.sm[i,],bandwidth = 5,x.points = 1:ncol(eeg50.sm))$y
}

# ===========================================================================
# THMM vs THSMM: Viterbi accuracy against the scored sleep stages
# ===========================================================================
#
# Run from anywhere:   Rscript osaModels/eegModels.R
# Environment variables: OSA_CORES (parallel workers, default detectCores() - 1)
# and OSA_QUICK = 1 (a small grid, to check that everything runs).
#
# Every combination of the following settings is one "setting":
#   nStates  2 (awake / asleep, eeg50.ss2), 3 (awake / REM / NREM, eeg50.ss3)
#            or 5 (the scored stages, eeg50.ss)
#   data     raw (eeg50.dat) or kernel-smoothed (eeg50.sm) log-power curves
#   band     frequency bins used: Freq 0-49 only (the curves are cut to these
#            50 bins before smoothing, so higher frequencies never enter)
#   norm     nonparametric mean-curve emissions with the L2, W21 or W22 norm,
#            log b_j(O) = -scale / 2 ||P(O - h_j)||^2
#   scale    "unit": scale = 1 (thmm1.R and tests/paperEeg.R); "auto": funHMM's
#            data-driven sigma (its default), converted to the package's scale.
#            The scale weighs the emissions against the dynamics, so it decides
#            how much a dwell model can change the Viterbi path.
#   start    "kmeans": k-means on the curves (unsupervised); "stages": the
#            scored-stage mean curves, transitions and mean dwells (an oracle
#            start: how well does each model hold on to the true structure?)
#
# Within a setting every model starts from the same mean curves, embedded
# transitions, mean dwell times and initial distribution, and evaluates the
# same emission functional, so their log-likelihoods and BICs are comparable
# within a setting (not across settings). Models:
#   THMM (funHMM)            funHMM::thmm(type = "nonpar"), sigma matched to the
#                            scale (checked: same log-likelihood and Viterbi path
#                            as the package's geometric THSMM at equal parameters)
#   THSMM, geometric         the package's THMM, a sanity check on funHMM
#   THSMM, Poisson / neg. binomial / Weibull     1 + Poisson, 1 + NB, discrete Weibull
#   THSMM, shifted Poisson   shifts profiled over shiftCandidates (2 states only)
#   THSMM, npg(m)            free dwell probabilities for 1, ..., m - 1 epochs
#                            plus a geometric tail, for several head lengths m
#
# Scores of each Viterbi path against the scored stages: accuracy and balanced
# accuracy (mean per-stage recall) after the best matching of fitted states to
# stages, and the ARI. The "winning" THSMM is the dwell family with the best
# mean accuracy rank across settings; the BIC table asks whether that winner
# could be picked without the stage labels.
#
# Output in osaModels/output/: results.csv (one row per fit), summary tables
# (.csv and LaTeX .tex), confusion.txt, figures (osaModels.pdf) and a cache of
# finished fits (cache/), so an interrupted run resumes where it stopped.
# Delete cache/ after changing the EM settings.

library(parallel)

nCores <- as.integer(Sys.getenv("OSA_CORES", max(1L, detectCores() - 1L)))
quick <- identical(Sys.getenv("OSA_QUICK"), "1")

outputDir <- file.path("osaModels", if (quick) "output-quick" else "output")
cacheDir <- file.path(outputDir, "cache")
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

truths <- list(`2` = eeg50.ss2, `3` = eeg50.ss3, `5` = eeg50.ss)
stageLabels <- list(`2` = c("Awake", "Asleep"), `3` = c("Awake", "REM", "NREM"),
                    `5` = c("Awake", "REM", "NREM1", "NREM2", "NREM3"))
datasets <- list(raw = eeg50.dat, smoothed = eeg50.sm)
bands <- list(`0-49` = 1:50)
norms <- c("L2", "W21", "W22")
scaleTypes <- c("unit", "auto")
startTypes <- c("kmeans", "stages")

parametricDwellSize <- 150        # exact dwell support; scored sojourns reach 169 epochs
shiftCandidates <- c(1, 2, 4, 8)
seed <- 137
emTol <- 1e-8
emMaxIter <- 1000

thmmMethod <- "THMM (funHMM)"
thsmmFamilies <- list(
  "THSMM, geometric"       = list(type = "geo"),
  "THSMM, Poisson"         = list(type = "pois", size = parametricDwellSize),
  "THSMM, neg. binomial"   = list(type = "nb", size = parametricDwellSize),
  "THSMM, Weibull"         = list(type = "weibull", size = parametricDwellSize),
  "THSMM, shifted Poisson" = list(type = "shiftp", size = parametricDwellSize, maxStates = 2),
  "THSMM, npg(10)"         = list(type = "npg", size = 10),
  "THSMM, npg(25)"         = list(type = "npg", size = 25),
  "THSMM, npg(50)"         = list(type = "npg", size = 50)
)

if (quick) {
  truths <- truths["3"]
  datasets <- datasets["smoothed"]
  norms <- "W21"
  startTypes <- "kmeans"
  thsmmFamilies <- thsmmFamilies[c("THSMM, geometric", "THSMM, neg. binomial", "THSMM, npg(10)")]
}

methodOrder <- c(thmmMethod, names(thsmmFamilies))
candidateFamilies <- setdiff(names(thsmmFamilies), "THSMM, geometric")   # competing for "winner"
settingColumns <- c("nStates", "data", "band", "norm", "scale", "start")


# ---------------------------------------------------------------------------
# Starting values and fits
# ---------------------------------------------------------------------------

# funHMM writes the functional as -||T(O - h)||^2 / (2 sigma^2) with T = 1 / sqrt(len)
# (L2), sqrt(len) * first differences (W21) or len^1.5 * second differences (W22);
# the package as -scale / 2 times the plain sum of squares.
normFactor <- function(norm, len) switch(norm, L2 = 1 / len, W21 = len, W22 = len^3)

emissionSettings <- function(x, norm, scaleType) {
  len <- ncol(x)
  sigma <- if (scaleType == "unit") {
    sqrt(normFactor(norm, len))
  } else {
    funHMM:::.sigma_auto(array(x, c(nrow(x), len, 1L)), "nonpar", norm, 0.5)
  }
  list(sigma = sigma, scale = normFactor(norm, len) / sigma^2)
}

# Shared start: mean curves, state labels, and the transitions / mean dwells /
# initial distribution implied by those labels.
makeStart <- function(x, nStates, control, startType, truth) {
  family <- getEmissionFamily("nonpar")
  prepared <- family$prepare(x, control)
  if (startType == "kmeans") {
    set.seed(seed)
    init <- family$initialize(prepared, nStates, control, "kmeans")
    means <- init$params
    labels <- init$labels
  } else {
    labels <- truth
    means <- family$update(prepared, labelsToWeights(labels, nStates), NULL, control,
                           list(minCount = 0, minSd = 1e-3))
  }
  delta <- tabulate(labels[1], nbins = nStates) + 0.1
  list(means = means, tpm = initializeEmbeddedTpmFromLabels(labels, nStates),
       meanDwells = initializeMeanDwellsFromLabels(labels, nStates), delta = delta / sum(delta))
}

fitFunHmm <- function(x, nStates, norm, emission, start) {
  stay <- 1 - pmin(pmax(1 / start$meanDwells, 1e-3), 1 - 1e-3)   # geometric dwell with the start's means
  trans <- diag(stay, nStates) + (1 - stay) * start$tpm
  fit <- suppressWarnings(funHMM::thmm(x, nStates, type = "nonpar", norm = norm, sigma = emission$sigma,
                                       par = start$means, trans = trans, init.prob = start$delta,
                                       tol = emTol, max.iter = emMaxIter))
  list(states = fit$states, logLik = fit$loglik, converged = fit$converged, iterations = fit$iter,
       nParameters = length(start$means) + nStates * (nStates - 1) + (nStates - 1),
       means = fit$par, trans = fit$trans)
}

fitThsmm <- function(x, nStates, control, start, family) {
  sojournType <- rep(family$type, nStates)
  sojournParameters <- lapply(seq_len(nStates), function(i) {
    initializeSojournOne(family$type, dwellSize = family$size %||% 1, meanDwell = start$meanDwells[i])
  })
  model <- buildHsmm(start$tpm, sojournParameters, sojournType, start$means, "nonpar",
                     delta = start$delta, stationary = FALSE, emissionControl = control)
  fit <- if (family$type == "shiftp") {
    suppressWarnings(fitHsmmShiftGrid(model, x, rep(list(shiftCandidates), nStates), stationary = FALSE,
                                      fitFunction = fitHsmmEM, maxCombinations = length(shiftCandidates)^nStates,
                                      tol = emTol, maxIter = emMaxIter))
  } else {
    suppressWarnings(fitHsmmEM(model, x, stationary = FALSE, tol = emTol, maxIter = emMaxIter))
  }
  fitted <- fit$fittedModel
  dwellCount <- sum(mapply(getSojournParamCount, fitted$sojournParameters, fitted$sojournType,
                           MoreArgs = list(fitShift = FALSE)))
  list(states = hsmmViterbiEM(fitted, x)$states, logLik = fit$logLik,
       converged = fit$converged %||% NA, iterations = fit$iterations %||% NA,
       nParameters = length(fitted$emissionParameters) + dwellCount + nStates * max(nStates - 2, 0) +
         (nStates - 1) + if (family$type == "shiftp") nStates else 0,   # shifts are profiled
       model = fitted[c("tpm", "sojournType", "sojournParameters", "emissionParameters", "delta")],
       bestShift = fit$bestShift)
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
permutations <- lapply(setNames(names(truths), names(truths)), function(k) allPermutations(as.integer(k)))

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
# Jobs
# ---------------------------------------------------------------------------

settings <- expand.grid(start = startTypes, scale = scaleTypes, norm = norms, band = names(bands),
                        data = names(datasets), nStates = as.integer(names(truths)),
                        stringsAsFactors = FALSE)[, settingColumns]
jobs <- merge(settings, data.frame(method = methodOrder, stringsAsFactors = FALSE))
keep <- vapply(seq_len(nrow(jobs)), function(j) {
  family <- thsmmFamilies[[jobs$method[j]]]
  is.null(family$maxStates) || jobs$nStates[j] <= family$maxStates
}, logical(1))
jobs <- jobs[keep, ]
# Longest first: shift grids, then many states and long curves.
jobs <- jobs[order(-(jobs$method == "THSMM, shifted Poisson"), -jobs$nStates, -nchar(jobs$band)), ]
rownames(jobs) <- NULL

jobKey <- function(job) {
  key <- paste(c(job[settingColumns], job$method, emTol, emMaxIter, parametricDwellSize,
                 paste(shiftCandidates, collapse = "-")), collapse = "_")
  gsub("[^A-Za-z0-9._-]", "", gsub("[ ,()]+", "-", key))
}

runJob <- function(job) {
  cacheFile <- file.path(cacheDir, paste0(jobKey(job), ".rds"))
  if (file.exists(cacheFile)) return(readRDS(cacheFile))
  nStates <- job$nStates
  x <- datasets[[job$data]][, bands[[job$band]], drop = FALSE]
  truth <- truths[[as.character(nStates)]]
  elapsed <- system.time(fit <- tryCatch({
    emission <- emissionSettings(x, job$norm, job$scale)
    control <- list(norm = job$norm, scale = emission$scale)
    start <- makeStart(x, nStates, control, job$start, truth)
    out <- if (job$method == thmmMethod) {
      fitFunHmm(x, nStates, job$norm, emission, start)
    } else {
      fitThsmm(x, nStates, control, start, thsmmFamilies[[job$method]])
    }
    c(out, list(sigma = emission$sigma, scaleValue = emission$scale))
  }, error = function(e) list(error = conditionMessage(e))))[["elapsed"]]
  fit$seconds <- elapsed
  saveRDS(fit, cacheFile)
  fit
}

cat(sprintf("osaModels: %d epochs, %d settings, %d fits on %d cores%s\n", nrow(eeg50.dat), nrow(settings),
            nrow(jobs), nCores, if (quick) " (quick run)" else ""))
started <- Sys.time()
fits <- mclapply(seq_len(nrow(jobs)), function(j) runJob(jobs[j, ]), mc.cores = nCores, mc.preschedule = FALSE)
cat(sprintf("Fitted in %.1f minutes\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

isFailed <- function(fit) !is.list(fit) || !is.null(fit$error) || is.null(fit$states) || !is.finite(fit$logLik)

results <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  fit <- fits[[j]]
  nStates <- jobs$nStates[j]
  failed <- isFailed(fit)
  score <- if (failed) c(accuracy = NA, balancedAccuracy = NA, ARI = NA) else
    scorePath(fit$states, truths[[as.character(nStates)]], nStates)
  data.frame(
    jobs[j, c(settingColumns, "method")],
    t(score),
    logLik = if (failed) NA else fit$logLik,
    nParameters = if (failed) NA else fit$nParameters,
    BIC = if (failed) NA else -2 * fit$logLik + fit$nParameters * log(nrow(eeg50.dat)),
    converged = if (failed) NA else fit$converged,
    iterations = if (failed) NA else fit$iterations,
    seconds = round(if (is.list(fit)) fit$seconds %||% NA else NA, 1),
    error = if (!is.list(fit)) as.character(fit) else if (is.null(fit$error)) "" else fit$error,
    row.names = NULL
  )
}))
results <- results[do.call(order, c(unname(as.list(results[settingColumns])),
                                     list(match(results$method, methodOrder)))), ]
results$setting <- do.call(paste, c(unname(as.list(results[settingColumns])), sep = " | "))

# Within each setting: rank by accuracy (1 = best) and gain over the THMM.
bySetting <- split(seq_len(nrow(results)), results$setting)
results$rank <- NA_real_
results$gainOverThmm <- NA_real_
for (rows in bySetting) {
  results$rank[rows] <- rank(-results$accuracy[rows], na.last = "keep", ties.method = "min")
  thmm <- results$accuracy[rows][results$method[rows] == thmmMethod]
  if (length(thmm) == 1L) results$gainOverThmm[rows] <- results$accuracy[rows] - thmm
}
rownames(results) <- NULL

failures <- results[results$error != "" | is.na(results$accuracy), ]
if (nrow(failures)) {
  cat(sprintf("%d of %d fits failed, e.g.:\n", nrow(failures), nrow(results)))
  print(head(failures[, c("setting", "method", "error")], 10), row.names = FALSE)
}
unconverged <- sum(!results$converged, na.rm = TRUE)
if (unconverged) cat(sprintf("%d fits reached emMaxIter = %d without converging\n", unconverged, emMaxIter))

geoCheck <- merge(results[results$method == thmmMethod, c("setting", "logLik")],
                  results[results$method == "THSMM, geometric", c("setting", "logLik")], by = "setting")
if (nrow(geoCheck)) {
  cat(sprintf("THMM (funHMM) vs geometric THSMM: largest log-likelihood difference %.3g\n",
              max(abs(geoCheck$logLik.x - geoCheck$logLik.y), na.rm = TRUE)))
}

saveRDS(list(jobs = jobs, fits = fits, results = results), file.path(outputDir, "fits.rds"))
write.csv(results, file.path(outputDir, "results.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# Summary tables
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

# Per method, over the settings of each (nStates, start): mean scores, mean rank,
# how often it is the best model, and how often / by how much it beats the THMM.
summarizeMethods <- function(rows) {
  groups <- split(rows, rows[, c("nStates", "start", "method")], drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    data.frame(g[1, c("nStates", "start", "method")],
               settings = nrow(g), failed = sum(is.na(g$accuracy)),
               accuracy = mean(g$accuracy, na.rm = TRUE),
               balancedAccuracy = mean(g$balancedAccuracy, na.rm = TRUE),
               ARI = mean(g$ARI, na.rm = TRUE),
               meanRank = mean(g$rank, na.rm = TRUE),
               best = sum(g$rank == 1, na.rm = TRUE),
               gainOverThmm = mean(g$gainOverThmm, na.rm = TRUE),
               beatsThmm = mean(g$gainOverThmm > 0, na.rm = TRUE),
               row.names = NULL)
  }))
  out[order(out$nStates, out$start, out$meanRank), ]
}
methodSummary <- summarizeMethods(results)

# Per setting: the THMM against the best THSMM (geometric excluded).
settingBest <- do.call(rbind, lapply(bySetting, function(rows) {
  r <- results[rows, ]
  candidates <- r[r$method %in% candidateFamilies & !is.na(r$accuracy), ]
  thmm <- r[r$method == thmmMethod, ]
  if (!nrow(candidates)) return(NULL)
  best <- candidates[which.max(candidates$accuracy), ]
  bic <- candidates[which.min(candidates$BIC), ]
  data.frame(r[1, settingColumns],
             thmmAccuracy = if (nrow(thmm)) thmm$accuracy else NA,
             bestThsmm = best$method, bestAccuracy = best$accuracy,
             gain = best$accuracy - (if (nrow(thmm)) thmm$accuracy else NA),
             bicThsmm = bic$method, bicAccuracy = bic$accuracy,
             bicGain = bic$accuracy - (if (nrow(thmm)) thmm$accuracy else NA),
             row.names = NULL)
}))
settingBest <- settingBest[do.call(order, unname(as.list(settingBest[settingColumns]))), ]
rownames(settingBest) <- NULL

# The best single configurations for each number of states.
topConfigurations <- do.call(rbind, lapply(split(results, results$nStates), function(r) {
  r <- r[!is.na(r$accuracy), ]
  head(r[order(-r$accuracy), c(settingColumns, "method", "accuracy", "balancedAccuracy", "ARI")], 10)
}))
rownames(topConfigurations) <- NULL

# The winner: the competing dwell family with the best mean accuracy rank,
# per number of states and over all settings with the unsupervised start.
winnerTable <- function(rows) {
  rows <- rows[rows$method %in% candidateFamilies, ]
  ranks <- tapply(rows$rank, rows$method, mean, na.rm = TRUE)
  gains <- tapply(rows$gainOverThmm, rows$method, mean, na.rm = TRUE)
  beats <- tapply(rows$gainOverThmm > 0, rows$method, mean, na.rm = TRUE)
  data.frame(method = names(ranks), meanRank = as.numeric(ranks), gainOverThmm = as.numeric(gains[names(ranks)]),
             beatsThmm = as.numeric(beats[names(ranks)]), row.names = NULL)[order(ranks), ]
}
unsupervised <- results[results$start == "kmeans", ]
if (!nrow(unsupervised)) unsupervised <- results
winners <- do.call(rbind, c(
  lapply(split(unsupervised, unsupervised$nStates), function(r) cbind(nStates = r$nStates[1], head(winnerTable(r), 3))),
  list(cbind(nStates = "all", head(winnerTable(unsupervised), 3)))
))
rownames(winners) <- NULL

tables <- list(
  methodSummary = list(table = methodSummary, caption = paste(
    "Sleep-stage EEG (CF050, C4-M1): Viterbi scores of each model, averaged over the data / band / norm /",
    "scale settings; meanRank and best rank models by accuracy within a setting, and gainOverThmm and",
    "beatsThmm compare with funHMM's THMM from the same start.")),
  settingBest = list(table = settingBest, caption = paste(
    "Per setting: accuracy of funHMM's THMM, of the most accurate THSMM, and of the THSMM with the",
    "smallest BIC (chosen without the stage labels).")),
  topConfigurations = list(table = topConfigurations, caption = paste(
    "The ten most accurate configurations for each number of states.")),
  winners = list(table = winners, caption = paste(
    "Best THSMM dwell families by mean accuracy rank (k-means start), for each number of states and overall."))
)
for (name in names(tables)) {
  write.csv(tables[[name]]$table, file.path(outputDir, paste0(name, ".csv")), row.names = FALSE)
  writeLatex(tables[[name]]$table, file.path(outputDir, paste0(name, ".tex")), tables[[name]]$caption)
}

printTable <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(format(table, digits = 3), row.names = FALSE)
}
printTable("Method summary (averaged over settings)", methodSummary)
printTable("Top configurations", topConfigurations)
printTable("Winning THSMM dwell families (k-means start)", winners)
cat(sprintf("\nBIC picks the most accurate THSMM in %d of %d settings; the BIC-chosen THSMM beats the THMM in %d.\n",
            sum(settingBest$bicThsmm == settingBest$bestThsmm), nrow(settingBest),
            sum(settingBest$bicGain > 0, na.rm = TRUE)))
overallWinner <- winners[winners$nStates == "all", ][1, ]
cat(sprintf("Winning THSMM over all k-means-start settings: %s (mean rank %.2f, beats the THMM in %.0f%% of settings, mean gain %+.3f)\n",
            overallWinner$method, overallWinner$meanRank, 100 * overallWinner$beatsThmm, overallWinner$gainOverThmm))

# Confusion tables of the best THSMM and the THMM in each setting.
sink(file.path(outputDir, "confusion.txt"))
for (s in seq_len(nrow(settingBest))) {
  setting <- settingBest[s, ]
  nStates <- setting$nStates
  truth <- truths[[as.character(nStates)]]
  for (method in unique(c(thmmMethod, setting$bestThsmm))) {
    j <- which(do.call(paste, unname(as.list(jobs[settingColumns]))) == do.call(paste, unname(as.list(setting[settingColumns]))) &
                 jobs$method == method)
    if (length(j) != 1L || isFailed(fits[[j]])) next
    cat(sprintf("\n%s | %s\n", paste(unlist(setting[settingColumns]), collapse = " | "), method))
    print(table(stage = factor(truth, seq_len(nStates), stageLabels[[as.character(nStates)]]),
                viterbi = factor(alignToTruth(fits[[j]]$states, truth, nStates), seq_len(nStates),
                                 stageLabels[[as.character(nStates)]])))
  }
}
sink()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

methodColours <- setNames(c("#000000", "#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#009E73",
                            "#7FBF7B", "#1B7837")[seq_along(methodOrder)], methodOrder)

pdf(file.path(outputDir, "osaModels.pdf"), width = 11, height = 7)

# Accuracy gain over the THMM, by dwell family, start and number of states.
for (startType in startTypes) {
  par(mfrow = c(1, length(truths)), mar = c(10, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (k in names(truths)) {
    r <- results[results$nStates == as.integer(k) & results$start == startType & results$method != thmmMethod, ]
    r$method <- factor(r$method, setdiff(methodOrder, thmmMethod))
    if (all(is.na(r$gainOverThmm))) { plot.new(); next }
    boxplot(gainOverThmm ~ method, data = r, las = 2, xlab = "", ylab = "accuracy - THMM accuracy",
            col = methodColours[levels(r$method)], main = sprintf("%s states", k), drop = FALSE)
    abline(h = 0, lty = 2)
  }
  mtext(sprintf("Viterbi accuracy gain over funHMM's THMM across settings (%s start)", startType),
        outer = TRUE, font = 2)
}

# Hypnograms: scored stages, the THMM, and the most accurate THSMM, in the
# setting where that THSMM is most accurate (k-means start).
for (k in names(truths)) {
  nStates <- as.integer(k)
  rows <- settingBest[settingBest$nStates == nStates & settingBest$start == startTypes[1], ]
  if (!nrow(rows) || all(is.na(rows$bestAccuracy))) next
  setting <- rows[which.max(rows$bestAccuracy), ]
  truth <- truths[[k]]
  shown <- unique(c(thmmMethod, setting$bestThsmm))
  par(mfrow = c(length(shown) + 1, 1), mar = c(1, 6, 1.5, 1), oma = c(3, 0, 2, 0))
  plot(truth, type = "s", yaxt = "n", xlab = "", ylab = "", main = "Scored stages", ylim = c(0.5, nStates + 0.5))
  axis(2, at = seq_len(nStates), labels = stageLabels[[k]], las = 1, cex.axis = 0.7)
  for (method in shown) {
    j <- which(do.call(paste, unname(as.list(jobs[settingColumns]))) == do.call(paste, unname(as.list(setting[settingColumns]))) &
                 jobs$method == method)
    if (length(j) != 1L || isFailed(fits[[j]])) { plot.new(); next }
    accuracy <- results$accuracy[results$setting == paste(unlist(setting[settingColumns]), collapse = " | ") &
                                   results$method == method]
    plot(alignToTruth(fits[[j]]$states, truth, nStates), type = "s", col = methodColours[[method]], yaxt = "n",
         xlab = "", ylab = "", ylim = c(0.5, nStates + 0.5), main = sprintf("%s (accuracy %.3f)", method, accuracy))
    axis(2, at = seq_len(nStates), labels = stageLabels[[k]], las = 1, cex.axis = 0.7)
  }
  mtext(sprintf("%d states: %s", nStates, paste(unlist(setting[settingColumns[-1]]), collapse = ", ")),
        outer = TRUE, font = 2)
  mtext("epoch", side = 1, outer = TRUE, line = 1.5)

  # Fitted dwell pmfs of the winning THSMM against the scored sojourn lengths.
  j <- which(do.call(paste, unname(as.list(jobs[settingColumns]))) == do.call(paste, unname(as.list(setting[settingColumns]))) &
               jobs$method == setting$bestThsmm)
  fit <- fits[[j]]
  if (length(j) == 1L && !isFailed(fit)) {
    runs <- rle(truth)
    aligned <- alignToTruth(fit$states, truth, nStates)
    stateOf <- vapply(seq_len(nStates), function(s) {
      matched <- unique(fit$states[aligned == s])
      if (length(matched)) matched[1] else NA_integer_
    }, integer(1))
    par(mfrow = c(1, nStates), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (s in seq_len(nStates)) {
      observed <- runs$lengths[runs$values == s]
      breaks <- seq(0, max(observed, 10) + 10, by = 5)
      hist(observed, breaks = breaks, freq = FALSE, col = "grey85", border = "white", xlab = "sojourn (epochs)",
           main = sprintf("%s (n = %d)", stageLabels[[k]][s], length(observed)))
      if (!is.na(stateOf[s])) {
        model <- fit$model
        pmf <- getModelDwellPmf(buildHsmm(model$tpm, model$sojournParameters, model$sojournType,
                                          model$emissionParameters, "nonpar", delta = model$delta,
                                          stationary = FALSE, emissionControl = list(norm = setting$norm)),
                                stateOf[s], max(breaks))
        lines(seq_along(pmf), pmf, col = methodColours[[setting$bestThsmm]], lwd = 2)
      }
    }
    mtext(sprintf("%d states: scored sojourns and the dwell pmf of %s", nStates, setting$bestThsmm),
          outer = TRUE, font = 2)
  }
}

invisible(dev.off())
cat(sprintf("\nResults written to %s\n", outputDir))
