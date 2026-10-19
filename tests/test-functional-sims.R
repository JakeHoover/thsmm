# Arbitrary dwell-time topological HMM: validation checks and an
# articleSims-style simulation study.
#
# Run from the package root:   Rscript tests/test-functional-sims.R
# (about 2-3 minutes). Figures and a results table are written to
# tests/output/.
#
# Part A -- correctness (PASS / FAIL; non-zero exit status on failure)
#   A1. With geometric dwell times, fitHsmmEM() reproduces thmm1.R's
#       Baum-Welch iterates exactly (nonparametric emissions, L2 / W21 / W22).
#       Needs thmm1.R saved as tests/thmm1.r; skipped otherwise.
#   A2. The closed-form emission M-steps (BMWD, fBM, 2-d BM, OU including the
#       rate >= 0 boundary, nonparametric L2 / W21 / W22) equal a numerical
#       maximization of the same weighted functional.
#   A3. The sparse E-step equals a dense log-space forward pass on
#       hsmmToHmm()'s transition matrix, including log-densities of order
#       -1e5, where the package's scaledForward() underflows.
#   A4. The sparse log-space Viterbi equals the package's dense viterbiLog().
#   A5. Every EM fit in Part B is monotone.
#
# Part B -- simulation study mirroring articleSims.r (Kashlak et al., 2023,
# Sec. 5): five states, 100-point curves on [0, 1], uniform embedded
# transitions. Data come from a semi-Markov truth (shifted-Poisson dwell
# times), except one control scenario with the paper's geometric A1 matrix.
# Each scenario is fitted with
#   * THMM, geometric dwell   -- the model class of thmm1.R (with the
#                                paper's functionals)
#   * HSMM-THMM, npg(10)      -- nonparametric dwell head + geometric tail
#   * HSMM-THMM, shiftp       -- true dwell family (shift fixed at the truth)
#   * fPCA(1) + Gaussian HMM  -- stand-in for the paper's fPCA + mhsmm
#                                baseline (mhsmm is not needed)
# and scored by ARI and cross-entropy of the Viterbi path (as in
# articleSims.r), BIC, dwell-pmf recovery, and parameter estimates.

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

outputDir <- file.path("tests", "output")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

nCurves <- 300                    # articleSims.r uses 200; more sojourns help dwell estimation
nPoints <- 100
grid <- seq(0, 1, length.out = nPoints)
K <- 5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

checks <- data.frame(test = character(), value = numeric(), limit = numeric(), pass = logical())

check <- function(test, value, limit) {
  pass <- is.finite(value) && value <= limit
  checks[nrow(checks) + 1L, ] <<- list(test, value, limit, pass)
  cat(sprintf("  [%s] %-62s %10.3e  (limit %.0e)\n", if (pass) "PASS" else "FAIL", test, value, limit))
}

# Cross-entropy of articleSims.r.
crossEnt <- function(est, tru) {
  tab <- table(est, tru)
  ftab <- t(t(tab) / colSums(tab))
  -sum(apply(log2(ftab), 2, max))
}

adjustedRand <- function(a, b) {
  if (requireNamespace("mclust", quietly = TRUE)) return(mclust::adjustedRandIndex(a, b))
  tab <- table(a, b)
  choose2 <- function(n) n * (n - 1) / 2
  index <- sum(choose2(tab))
  rows <- sum(choose2(rowSums(tab))); cols <- sum(choose2(colSums(tab)))
  expected <- rows * cols / choose2(sum(tab))
  (index - expected) / (0.5 * (rows + cols) - expected)
}

permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  smaller <- permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))))
}
allPermutations <- permutations(K)

# Estimated state that best matches each true state (maximizes matched counts).
matchStates <- function(estimated, truth) {
  tab <- table(factor(estimated, seq_len(K)), factor(truth, seq_len(K)))
  scores <- apply(allPermutations, 1, function(p) sum(tab[cbind(p, seq_len(K))]))
  allPermutations[which.max(scores), ]
}

countParameters <- function(model) {
  dwell <- sum(mapply(getSojournParamCount, model$sojournParameters, model$sojournType,
                      MoreArgs = list(fitShift = FALSE)))
  length(model$emissionParameters) + dwell + max(K * (K - 2), 0) + (K - 1)
}

maxDwell <- 40

# Fit, decode, and score one model against the truth.
fitAndScore <- function(method, start, x, truth, trueStates, scenario) {
  elapsed <- system.time(
    fit <- suppressWarnings(fitHsmmEM(start, x, stationary = FALSE, tol = 1e-9, maxIter = 3000))
  )[["elapsed"]]
  decoded <- hsmmViterbiEM(fit$fittedModel, x)
  match <- matchStates(decoded$states, trueStates)
  model <- fit$fittedModel

  dwellError <- mean(vapply(seq_len(K), function(i) {
    sum(abs(getModelDwellPmf(model, match[i], maxDwell) - getModelDwellPmf(truth, i, maxDwell)))
  }, numeric(1)))

  exact <- getEmissionFamily(model$emissionType)$exactMStep(model$emissionControl)
  if (exact) {
    check(sprintf("%s / %s: EM monotone", scenario, method), max(0, -diff(fit$logLikTrace)), 1e-7)
  }

  list(
    row = data.frame(
      scenario = scenario, method = method,
      ARI = adjustedRand(decoded$states, trueStates),
      Ent = crossEnt(decoded$states, trueStates),
      logLik = fit$logLik, nPar = countParameters(model),
      BIC = -2 * fit$logLik + countParameters(model) * log(nrow(as.matrix(x))),
      dwellL1 = dwellError, iterations = fit$iterations, seconds = elapsed
    ),
    fit = fit, decoded = decoded, match = match
  )
}

# Same emission start for every dwell family: initializeHsmmEM() is
# deterministic given the seed, so only the dwell starting values differ.
startModel <- function(x, sojournType, emissionType, control, truth = NULL, seed = 137) {
  init <- initializeHsmmEM(x, K, sojournType, emissionType, control, dwellSize = 10,
                           stationary = FALSE, method = "kmeans", seed = seed)

  if (sojournType == "shiftp" && !is.null(truth)) {
    # Oracle shift: fixed at the truth (fitHsmmShiftGrid would profile it). The
    # initializer orders states by emission parameter, so true state i is first
    # matched to its starting state by parameter distance.
    model <- init$model
    costs <- apply(allPermutations, 1, function(p) {
      sum((model$emissionParameters[p, , drop = FALSE] - truth$emissionParameters)^2)
    })
    toStart <- allPermutations[which.min(costs), ]
    for (i in seq_len(K)) {
      j <- toStart[i]
      shift <- truth$sojournParameters[[i]][["shift"]]
      model$sojournParameters[[j]] <- c(lambda = max(init$meanDwells[j] - shift, 0.5),
                                        shift = shift, dwellSize = 20)
    }
    return(do.call(buildHsmm, model[c("tpm", "sojournParameters", "sojournType", "emissionParameters",
                                      "emissionType", "delta", "stationary", "fitShift", "emissionControl")]))
  }

  init$model
}

# fPCA(1) + Gaussian HMM baseline.
fpcaBaseline <- function(x, trueStates, scenario) {
  flat <- matrix(x, nrow = dim(x)[1])
  scores <- prcomp(flat)$x[, 1]
  init <- initializeHsmmEM(scores, K, "geo", "normal", dwellSize = 1, stationary = FALSE, seed = 137)
  fit <- suppressWarnings(fitHsmmEM(init$model, scores, stationary = FALSE, tol = 1e-9))
  decoded <- hsmmViterbiEM(fit$fittedModel, scores)
  data.frame(scenario = scenario, method = "fPCA(1) + Gaussian HMM",
             ARI = adjustedRand(decoded$states, trueStates), Ent = crossEnt(decoded$states, trueStates),
             logLik = NA, nPar = NA, BIC = NA, dwellL1 = NA, iterations = fit$iterations, seconds = NA)
}

omega <- matrix(1 / (K - 1), K, K); diag(omega) <- 0

# Semi-Markov truth: shifted-Poisson dwell times, means 3, 5, 5, 6, 4.
semiMarkovDwell <- Map(function(lambda, shift) c(lambda = lambda, shift = shift, dwellSize = 20),
                       c(1, 2, 1, 3, 2), c(2, 3, 4, 3, 2))

makeTruth <- function(emissionType, params, control, geometric = FALSE) {
  if (geometric) {
    # The paper's A1: stay 0.64, move 0.09 each  <=>  geo(0.36) dwell, uniform omega.
    sojourn <- rep(list(c(prob = 0.36)), K); type <- rep("geo", K)
  } else {
    sojourn <- semiMarkovDwell; type <- rep("shiftp", K)
  }
  buildHsmm(omega, sojourn, type, params, emissionType, delta = c(1, rep(0, K - 1)),
            stationary = FALSE, emissionControl = control)
}

results <- list()
figureData <- list()

runScenario <- function(scenario, truth, fitControl = truth$emissionControl, fits = c("geo", "npg", "shiftp"),
                        label = scenario, seed = 137, baseline = TRUE) {
  sim <- simulateHsmm(truth, nCurves, initialType = "aggregateStart", seed = seed)
  out <- list()

  for (dwell in fits) {
    method <- switch(dwell, geo = "THMM, geometric dwell", npg = "HSMM-THMM, npg(10)",
                     shiftp = "HSMM-THMM, shiftp (oracle shift)")
    start <- startModel(sim$x, dwell, truth$emissionType, fitControl, truth)
    out[[dwell]] <- fitAndScore(method, start, sim$x, truth, sim$states, label)
  }

  table <- do.call(rbind, lapply(out, `[[`, "row"))
  if (baseline) table <- rbind(table, fpcaBaseline(sim$x, sim$states, label))
  results[[label]] <<- table
  figureData[[label]] <<- list(sim = sim, truth = truth, out = out)
  invisible(list(sim = sim, out = out, table = table))
}

printScenario <- function(label, parameterLabels, transformParams = identity) {
  data <- figureData[[label]]
  table <- results[[label]]
  cat("\n------------------------------------------------------------------\n")
  cat(label, "\n")
  print(format(table[, c("method", "ARI", "Ent", "logLik", "nPar", "BIC", "dwellL1", "iterations")],
               digits = 4), row.names = FALSE)

  best <- data$out[["npg"]]
  if (is.null(best)) return(invisible())
  cat("\n  HSMM-THMM npg(10): Viterbi states (rows) vs true states (columns)\n")
  print(table(estimated = factor(best$decoded$states, best$match, letters[1:K]),
              truth = factor(data$sim$states, 1:K, LETTERS[1:K])))

  estimated <- transformParams(best$fit$fittedModel$emissionParameters)[best$match, , drop = FALSE]
  truthParams <- transformParams(data$truth$emissionParameters)
  meanDwell <- function(model, states) vapply(states, function(i) sum(seq_len(maxDwell) * getModelDwellPmf(model, i, maxDwell)), 1)
  comparison <- rbind(t(truthParams), t(estimated),
                      meanDwell(data$truth, 1:K), meanDwell(best$fit$fittedModel, best$match))
  rownames(comparison) <- c(paste("true", parameterLabels), paste("est ", parameterLabels),
                            "true mean dwell", "est  mean dwell")
  colnames(comparison) <- LETTERS[1:K]
  print(round(comparison, 2))
}


# ===========================================================================
# Part A: correctness
# ===========================================================================

cat("\nA1. Geometric dwell reproduces thmm1.R's Baum-Welch exactly\n")

thmmPath <- file.path("tests", "thmm1.r")

if (file.exists(thmmPath)) {
  suppressMessages(source(thmmPath, local = (thmmEnv <- new.env())))
  with(thmmEnv, {
    set.seed(137)
    tmat1 <- matrix(0.09, 5, 5) + 0.55 * diag(5)
    tt <- seq(0, 1, length.out = 100)
    dat <- sim.thmm(n = 200, prInit = c(1, 0, 0, 0, 0), prTran = tmat1, type = "nonpar", len = 100,
                    drift = t(sapply(0:4, function(j) sin(2 * pi * (tt + 0.2 * j)))))
    startMeans <- dat$data[sample(200, 5), ]
    A0 <- rdirichlet(5, rep(6, 5)); eta0 <- as.vector(rdirichlet(1, rep(6, 5)))
  })
  nIter <- 40

  for (norm in c("L2", "W21", "W22")) {
    th <- with(thmmEnv, thmm.init(dat$data, startMeans, A0, eta0, type = "nonpar"))
    th$prInit <- log(th$prInit); th$prTran <- log(th$prTran)
    like <- numeric(nIter)
    for (it in seq_len(nIter)) {
      th <- thmmEnv$thmm.bj.all(th, log.scale = TRUE, type = "nonpar", norm = norm)
      th <- thmmEnv$thmm.alfa(th, TRUE, "nonpar")
      th <- thmmEnv$thmm.beta(th, TRUE, "nonpar")
      like[it] <- thmmEnv$sumLogs.vec(th$alfa[nrow(th$alfa), ])
      th <- thmmEnv$thmm.rest(th, TRUE, type = "nonpar")
    }

    A0 <- thmmEnv$A0
    p <- 1 - diag(A0); omega0 <- A0 / p; diag(omega0) <- 0
    m0 <- buildHsmm(omega0, lapply(p, function(q) c(prob = q)), rep("geo", 5), thmmEnv$startMeans, "nonpar",
                    delta = thmmEnv$eta0, stationary = FALSE, emissionControl = list(norm = norm))
    em <- suppressWarnings(fitHsmmEM(m0, thmmEnv$dat$data, stationary = FALSE, maxIter = nIter, tol = -1))
    hazard <- unlist(em$fittedModel$sojournParameters)
    Aem <- diag(1 - hazard) + em$fittedModel$tpm * hazard

    check(sprintf("nonpar %s: max |log L trace, EM - thmm1| (40 its)", norm),
          max(abs(like - em$logLikTrace[seq_len(nIter)])), 1e-8)
    check(sprintf("nonpar %s: max |A, EM - thmm1|", norm), max(abs(exp(th$prTran) - Aem)), 1e-10)
    check(sprintf("nonpar %s: max |mean curves, EM - thmm1|", norm),
          max(abs(th$stPar - em$fittedModel$emissionParameters)), 1e-10)
  }
} else {
  cat("  (skipped: save thmm1.R as tests/thmm1.r to run this check)\n")
}


cat("\nA2. Closed-form emission M-steps maximize the weighted functional\n")

set.seed(11)
stateWeights <- matrix(rexp(nCurves * K), nCurves, K)
stateWeights <- stateWeights / rowSums(stateWeights)

weightedFunctional <- function(family, prepared, params, control, j) {
  sum(stateWeights[, j] * family$logDensity(prepared, params, control)[, j])
}

mStepGap <- function(label, type, params, control, lower = -Inf, nGrid = nPoints) {
  truth <- makeTruth(type, params, c(control, list(grid = seq(0, 1, length.out = nGrid))))
  x <- simulateHsmm(truth, nCurves, initialType = "aggregateStart", seed = 3)$x
  family <- getEmissionFamily(type)
  prepared <- family$prepare(x, truth$emissionControl)
  closedForm <- family$update(prepared, stateWeights, params, truth$emissionControl, list(minCount = 0, minSd = 0))

  gap <- 0
  for (j in seq_len(K)) {
    objective <- function(theta) {
      candidate <- closedForm; candidate[j, ] <- theta
      -weightedFunctional(family, prepared, candidate, truth$emissionControl, j)
    }
    numeric <- optim(closedForm[j, ] + 0.3, objective, method = "L-BFGS-B", lower = lower,
                     control = list(factr = 1, pgtol = 1e-12, maxit = 5000))
    # Positive gap = optim found a better point than the closed form. (Objective
    # values, not parameters, are compared: W21 / W22 are seminorms, so their
    # maximizing mean curve is only defined up to a constant / linear term.)
    gap <- max(gap, objective(closedForm[j, ]) - numeric$value)
  }
  check(sprintf("%s: closed form vs optim (objective gap)", label), gap, 1e-5)
}

mStepGap("BMWD", "bmwd", matrix(c(-4, -2, 0, 2, 4)), list())
mStepGap("fBM, hurst 0.8", "bmwd", matrix(c(-10, -6, -2, 0, 2)), list(hurst = 0.8))
mStepGap("fBM, hurst 0.3", "bmwd", matrix(c(-10, -6, -2, 0, 2)), list(hurst = 0.3))
mStepGap("2-d BM with drift", "bmwd", 3 * cbind(cos(2 * pi * (1:K) / K), sin(2 * pi * (1:K) / K)), list())
mStepGap("OU", "ou", cbind(c(4, 4, 8, 2, 20) * c(-2, 0, 4, 2, 1), c(4, 4, 8, 2, 20)), list(),
         lower = c(-Inf, 0))
# Explosive paths (negative true rate): the M-step must stop at rate = 0.
mStepGap("OU at the rate >= 0 boundary", "ou", cbind(rep(1, K), rep(-3, K)), list(), lower = c(-Inf, 0))
for (norm in c("L2", "W21", "W22")) {
  mStepGap(paste("nonpar", norm), "nonpar",
           t(sapply(0:4, function(j) sin(2 * pi * (seq(0, 1, length.out = 15) + 0.2 * j)))),
           list(norm = norm), nGrid = 15)
}


cat("\nA3. Sparse E-step vs dense log-space forward pass (large |log-density|)\n")

denseLogForward <- function(model, x) {
  hmm <- hsmmToHmm(model, initialType = "aggregateStart")
  family <- getEmissionFamily(model$emissionType)
  logB <- family$logDensity(family$prepare(x, model$emissionControl), model$emissionParameters,
                            model$emissionControl)[, rep(seq_len(K), hmm$aggregateSizes)]
  a <- log(hmm$delta) + logB[1, ]
  for (t in 2:nrow(logB)) {
    m <- max(a)
    a <- log(as.vector(exp(a - m) %*% hmm$gamma)) + m + logB[t, ]
  }
  m <- max(a)
  list(logLik = m + log(sum(exp(a - m))), logB = logB, hmm = hmm)
}

tt <- grid
sharpTruth <- makeTruth("nonpar", t(sapply(0:4, function(j) sin(2 * pi * (tt + 0.2 * j)))),
                        list(grid = grid, norm = "L2", scale = 2000))
sharpX <- simulateHsmm(sharpTruth, nCurves, initialType = "aggregateStart", seed = 5)$x
sharpStart <- startModel(sharpX, "npg", "nonpar", sharpTruth$emissionControl)
dense <- denseLogForward(sharpStart, sharpX)
sparse <- hsmmEStep(sharpStart, sharpX, stationary = FALSE)
naive <- scaledForward(dense$hmm$delta, dense$hmm$gamma, pmax(exp(dense$logB), 0))

cat(sprintf("  log-densities range from %.0f to %.0f; scaledForward() logLik = %s\n",
            min(dense$logB), max(dense$logB), format(naive$logLik)))
check("nonpar, scale 2000: |sparse - dense| logLik", abs(sparse$logLik - dense$logLik), 1e-6)

for (case in list(list("bmwd", matrix(c(-4, -2, 0, 2, 4)), list(grid = grid)),
                  list("ou", cbind(c(4, 4, 8, 2, 20) * c(-2, 0, 4, 2, 1), c(4, 4, 8, 2, 20)), list(grid = grid)))) {
  truth <- makeTruth(case[[1]], case[[2]], case[[3]])
  x <- simulateHsmm(truth, nCurves, initialType = "aggregateStart", seed = 5)$x
  start <- startModel(x, "npg", case[[1]], case[[3]])
  check(sprintf("%s: |sparse - dense| logLik", case[[1]]),
        abs(hsmmEStep(start, x, FALSE)$logLik - denseLogForward(start, x)$logLik), 1e-6)
}


cat("\nA4. Sparse log-space Viterbi vs the package's dense viterbiLog()\n")

viterbiCheck <- function(label, model, x) {
  hmm <- hsmmToHmm(model, initialType = "aggregateStart")
  family <- getEmissionFamily(model$emissionType)
  logB <- family$logDensity(family$prepare(x, model$emissionControl), model$emissionParameters,
                            model$emissionControl)[, rep(seq_len(K), hmm$aggregateSizes)]
  # Viterbi paths are invariant to per-time constants, so dense viterbiLog() can be
  # fed offset densities exp(logB - rowMax) and its log-probability corrected.
  offset <- apply(logB, 1, max)
  denseResult <- viterbiLog(hmm$delta, hmm$gamma, exp(logB - offset))
  sparseResult <- hsmmViterbiEM(model, x, stationary = FALSE)
  check(paste(label, ": expanded paths differ at # times"),
        sum(denseResult$expandedStates != sparseResult$expandedStates), 0)
  check(paste(label, ": |log-probability difference|"),
        abs(denseResult$logProb + sum(offset) - sparseResult$logProb), 1e-6)
}

viterbiCheck("nonpar npg", sharpStart, sharpX)
bmTruth <- makeTruth("bmwd", matrix(c(-4, -2, 0, 2, 4)), list(grid = grid))
bmX <- simulateHsmm(bmTruth, nCurves, initialType = "aggregateStart", seed = 8)$x
viterbiCheck("bmwd shiftp (truth)", bmTruth, bmX)
normalModel <- buildHsmm(omega, semiMarkovDwell, rep("shiftp", K), cbind(c(-4, -2, 0, 2, 4), 1), "normal",
                         delta = c(1, rep(0, K - 1)), stationary = FALSE)
viterbiCheck("normal shiftp", normalModel, simulateHsmm(normalModel, 400, "aggregateStart", seed = 2)$x)


# ===========================================================================
# Part B: simulation study (articleSims.r analogue)
# ===========================================================================

cat("\n\nB. Simulation study (A5 monotonicity checks are reported as fits run)\n")

# 5.1 Brownian motion with drift: low and medium separation
runScenario("BMWD, low separation", makeTruth("bmwd", matrix(c(-4, -2, 0, 2, 4)), list(grid = grid)))
runScenario("BMWD, medium separation", makeTruth("bmwd", matrix(c(-8, -4, 0, 4, 8)), list(grid = grid)))
# Control: the paper's geometric A1 truth (oracle = geometric, so no shiftp fit)
runScenario("BMWD, low separation, geometric truth (A1)",
            makeTruth("bmwd", matrix(c(-4, -2, 0, 2, 4)), list(grid = grid), geometric = TRUE),
            fits = c("geo", "npg"))

# 5.1 (thmm1's BM2d): 2-d Brownian motion, drift vectors on a circle
runScenario("2-d BM with drift",
            makeTruth("bmwd", 3 * cbind(cos(2 * pi * (1:K) / K), sin(2 * pi * (1:K) / K)), list(grid = grid)))

# 5.2 Ornstein-Uhlenbeck: mu = (-2, 0, 4, 2, 1), c = (4, 4, 8, 2, 20); b0 = c mu, b1 = c
ouRate <- c(4, 4, 8, 2, 20)
runScenario("OU process", makeTruth("ou", cbind(ouRate * c(-2, 0, 4, 2, 1), ouRate), list(grid = grid)))

# 5.3 Fractional BM, data hurst 0.8 and 0.25, fitted with hurst 0.25 / 0.5 / 0.8 (Tables 5-6)
for (setting in list(list(hurst = 0.8, drift = c(-10, -6, -2, 0, 2)),
                     list(hurst = 0.25, drift = c(-16, -8, -4, 0, 8)))) {
  truth <- makeTruth("bmwd", matrix(setting$drift), list(grid = grid, hurst = setting$hurst))
  for (fitHurst in c(0.25, 0.5, 0.8)) {
    runScenario("fBM", truth, fitControl = list(grid = grid, hurst = fitHurst), fits = c("geo", "npg"),
                label = sprintf("fBM, data hurst %.2f, fitted hurst %.2f", setting$hurst, fitHurst),
                baseline = fitHurst == 0.5)
  }
}

# 5.4 Nonparametric: phase-shifted sinusoids + Brownian-bridge K-L noise, L2 and W21 norms
sinusoids <- t(sapply(0:4, function(j) sin(2 * pi * (grid + 0.2 * j))))
for (norm in c("L2", "W21")) {
  runScenario(paste("Nonparametric sinusoids,", norm, "norm"),
              makeTruth("nonpar", sinusoids, list(grid = grid, norm = norm)),
              fits = c("geo", "npg", "shiftp"))
}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

identityParams <- function(p) p
printScenario("BMWD, low separation", "drift")
printScenario("BMWD, medium separation", "drift")
printScenario("BMWD, low separation, geometric truth (A1)", "drift")
printScenario("2-d BM with drift", c("drift x", "drift y"))
printScenario("OU process", c("mean", "rate"),
              transformParams = function(p) cbind(mean = ouMeanLevel(p)[, 1], rate = p[, 2]))
for (label in grep("^fBM", names(results), value = TRUE)) {
  cat("\n------------------------------------------------------------------\n", label, "\n")
  data <- figureData[[label]]
  for (dwell in names(data$out)) {
    fitDrift <- data$out[[dwell]]$fit$fittedModel$emissionParameters[data$out[[dwell]]$match, 1]
    entry <- data$out[[dwell]]$row
    cat(sprintf("  %-24s ARI %.3f  Ent %.3f  drift %s\n", entry$method, entry$ARI, entry$Ent,
                paste(sprintf("%7.2f", fitDrift), collapse = "")))
  }
  cat(sprintf("  %-46s %s  (truth)\n", "", paste(sprintf("%7.2f", data$truth$emissionParameters[, 1]), collapse = "")))
}
for (label in grep("^Nonparametric", names(results), value = TRUE)) {
  cat("\n------------------------------------------------------------------\n", label, "\n")
  print(format(results[[label]][, c("method", "ARI", "Ent", "logLik", "dwellL1", "iterations")], digits = 4),
        row.names = FALSE)
}

summary <- do.call(rbind, results)
rownames(summary) <- NULL
write.csv(summary, file.path(outputDir, "functional-sims-results.csv"), row.names = FALSE)

cat("\n\n==================================================================\n")
cat("Summary: ARI of the Viterbi path by scenario and method\n")
wide <- reshape(summary[, c("scenario", "method", "ARI")], idvar = "scenario", timevar = "method", direction = "wide")
names(wide) <- sub("^ARI\\.", "", names(wide))
print(format(wide, digits = 3), row.names = FALSE)


# ---------------------------------------------------------------------------
# Figures: decoded curves and dwell-time recovery
# ---------------------------------------------------------------------------

pdf(file.path(outputDir, "functional-sims.pdf"), width = 11, height = 6.5)
palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

for (label in names(figureData)) {
  data <- figureData[[label]]
  if (is.null(data$out$npg)) next
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

  firstCoordinate <- if (length(dim(data$sim$x)) == 3L) data$sim$x[, , 1] else data$sim$x
  matplot(grid, t(firstCoordinate), type = "l", lty = 1, col = palette[data$sim$states],
          xlab = expression(tau), ylab = "O(tau)", main = "Curves, coloured by true state")

  for (i in seq_len(K)) {
    truePmf <- getModelDwellPmf(data$truth, i, 15)
    fitted <- lapply(intersect(c("geo", "npg", "shiftp"), names(data$out)), function(dwell) {
      entry <- data$out[[dwell]]
      list(dwell = dwell, pmf = getModelDwellPmf(entry$fit$fittedModel, entry$match[i], 15))
    })
    upper <- max(truePmf, unlist(lapply(fitted, `[[`, "pmf")))
    plot(seq_len(15), truePmf, type = "h", lwd = 6, col = "grey70",
         ylim = c(0, upper * 1.05), xlab = "dwell time", ylab = "probability",
         main = paste("Dwell distribution, state", LETTERS[i]))
    for (f in fitted) {
      lines(seq_len(15), f$pmf, type = "b", pch = 19, cex = 0.6,
            col = c(geo = "#2166ac", npg = "#b2182b", shiftp = "#1a9850")[f$dwell])
    }
    if (i == 1) legend("topright", c("truth", "geometric", "npg(10)", "shiftp"), bty = "n", cex = 0.8,
                       col = c("grey70", "#2166ac", "#b2182b", "#1a9850"), lwd = c(6, 1, 1, 1))
  }
  mtext(label, outer = TRUE, font = 2)
}
invisible(dev.off())

cat(sprintf("\nWrote %s and %s\n", file.path(outputDir, "functional-sims-results.csv"),
            file.path(outputDir, "functional-sims.pdf")))


# ---------------------------------------------------------------------------

cat(sprintf("\n%d / %d checks passed\n", sum(checks$pass), nrow(checks)))
if (!all(checks$pass)) quit(status = 1)
