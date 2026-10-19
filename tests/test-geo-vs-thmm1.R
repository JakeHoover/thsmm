# THSMM with geometric dwell times vs the THMM of thmm1.r, on simulated data.
#
# Run from the package root:   Rscript tests/test-geo-vs-thmm1.R
#
# An HSMM whose dwell times are all geometric is an HMM: with dwell
# d_i ~ geo(p_i) on {1, 2, ...} and embedded transition matrix omega
# (zero diagonal),
#
#   A = diag(1 - p) + diag(p) omega,   i.e.   p_i = 1 - A_ii,
#                                             omega_ij = A_ij / (1 - A_ii).
#
# So fitHsmmEM() with sojournType = "geo" should reproduce thmm1.r's
# Baum-Welch exactly, whenever both codes evaluate the same emission
# functional. thmm1.r's nonparametric emissions are the package's "nonpar"
# family; for BMWD / BM2d / OU / fBM, thmm1.r's E-step functional differs
# from its own M-step (see the header of emissionFamilies.R), so the exact
# checks below run thmm1.r's forward pass, backward pass, re-estimation
# (thmm.rest, including its drift M-step) and Viterbi with the emission
# matrix thmm$bj supplied by the package's family. Part 4 then fits each code
# with its own functional and compares them statistically.
#
#   1. Structure: a geo-dwell HSMM is the HMM with matrix A (hsmmToHmm,
#      parametrization round trip, simulated chains of both codes).
#   2. Fixed parameters: log-likelihood, state posteriors, and Viterbi path /
#      log-probability agree with thmm1.r's alfa / beta / viterbi.
#   3. Baum-Welch iterates from a common start (40 iterations): log-likelihood
#      trace, A, initial distribution and emission parameters agree.
#      nonpar (L2, W21, W22) and nonparMulti use thmm1.r unmodified; BMWD and
#      BM2d use thmm1.r with the package emission matrix. OU and fBM: one
#      iteration, where only the dynamics are comparable exactly (thmm1.r's OU
#      M-step maximizes another functional by optim(); its fBM M-step
#      discretizes the drift kernel differently, compared to a tolerance).
#   4. Full fits to convergence on data simulated by thmm1.r's sim.thmm():
#      nonpar L2 with thmm.BaumWelch() itself, then a small Monte Carlo study
#      for BMWD where each code uses its own functional: both recover the
#      truth, their decodings agree, and the package is not worse.

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

thmmPath <- file.path("tests", "thmm1.r")
if (!file.exists(thmmPath)) stop("thmm1.r not found at ", thmmPath)
thmm1 <- new.env()
suppressMessages(source(thmmPath, local = thmm1))   # also attaches mclust

K <- 5
nCurves <- 200
nPoints <- 100
nIter <- 40


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

checks <- data.frame(test = character(), value = numeric(), limit = numeric(), pass = logical())

check <- function(test, value, limit) {
  pass <- is.finite(value) && value <= limit
  checks[nrow(checks) + 1L, ] <<- list(test, value, limit, pass)
  cat(sprintf("  [%s] %-64s %10.3e  (limit %.2g)\n", if (pass) "PASS" else "FAIL", test, value, limit))
}

logSumExp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

# THMM (A, eta) <-> geo-dwell HSMM (omega, p, delta).
thmmToHsmm <- function(A, eta, emissionParameters, emissionType, emissionControl = NULL) {
  p <- 1 - diag(A)
  omega <- A / p
  diag(omega) <- 0
  buildHsmm(omega, lapply(p, function(q) c(prob = q)), rep("geo", nrow(A)), emissionParameters,
            emissionType, delta = as.vector(eta), stationary = FALSE, emissionControl = emissionControl)
}

hsmmToThmmA <- function(model) {
  hazard <- vapply(model$sojournParameters, function(s) s[["prob"]], numeric(1))
  diag(1 - hazard) + hazard * model$tpm
}

# thmm1.r stores parameters per type in its own layout; these convert them to
# the package's nStates x p matrix.
packageParams <- function(stPar, type) {
  switch(type,
         BMWD = matrix(stPar, ncol = 1),
         BM2d = t(stPar),
         OU = cbind(stPar[1, ], stPar[2, ]),
         nonpar = stPar,
         nonparMulti = matrix(stPar, nrow = dim(stPar)[1]))   # c(curve 1, curve 2, ...)
}

packageData <- function(data, type) if (type == "BM2d") aperm(data, c(2, 3, 1)) else data

packageType <- function(type) switch(type, BMWD = "bmwd", BM2d = "bmwd", OU = "ou", "nonpar")

# Package log-density matrix at thmm1.r's current parameters.
packageLogDensity <- function(th, type, control) {
  family <- getEmissionFamily(packageType(type))
  prepared <- family$prepare(packageData(th$data, type), control)
  family$logDensity(prepared, packageParams(th$stPar, type), control)
}

# thmm$bj from thmm1.r's own functional or from the package family.
setEmission <- function(th, type, norm, packageEmission, control) {
  if (!packageEmission) return(thmm1$thmm.bj.all(th, log.scale = TRUE, type = type, norm = norm))
  th$bj <- packageLogDensity(th, type, control)
  th
}

# thmm1.r's Baum-Welch loop (as in thmm.BaumWelch, without printing). With
# packageEmission = TRUE the emission matrix comes from the package family;
# everything else is thmm1.r's code.
runThmm1 <- function(data, stPar, A, eta, type, iterations, norm = "W21", hurst = 0.5,
                     packageEmission = FALSE, control = NULL) {
  th <- thmm1$thmm.init(data, stPar, A, eta, type = type)
  th$prInit <- log(th$prInit)
  th$prTran <- log(th$prTran)
  like <- numeric(iterations)
  for (it in seq_len(iterations)) {
    th <- setEmission(th, type, norm, packageEmission, control)
    th <- thmm1$thmm.alfa(th, TRUE, type)
    th <- thmm1$thmm.beta(th, TRUE, type)
    like[it] <- thmm1$sumLogs.vec(th$alfa[nrow(th$alfa), ])
    th <- thmm1$thmm.rest(th, TRUE, type = type, hurst = hurst)
  }
  th$like <- like
  th
}

# Run thmm1.r and the package from the same start and compare the iterates.
compareIterates <- function(label, data, stPar, A, eta, type, iterations = nIter, norm = "W21",
                            hurst = 0.5, packageEmission = FALSE, control = NULL,
                            compareEmission = TRUE, emissionLimit = 1e-10) {
  th <- runThmm1(data, stPar, A, eta, type, iterations, norm, hurst, packageEmission, control)
  start <- thmmToHsmm(A, eta, packageParams(stPar, type), packageType(type), control)
  em <- suppressWarnings(fitHsmmEM(start, packageData(data, type), stationary = FALSE,
                                   maxIter = iterations, tol = -1))
  fitted <- em$fittedModel

  check(sprintf("%s: max |log L trace, EM - thmm1| (%d its)", label, iterations),
        max(abs(th$like - em$logLikTrace[seq_len(iterations)])), 1e-8)
  check(sprintf("%s: max |A, EM - thmm1|", label), max(abs(exp(th$prTran) - hsmmToThmmA(fitted))), 1e-10)
  check(sprintf("%s: max |initial dist, EM - thmm1|", label), max(abs(exp(th$prInit) - fitted$delta)), 1e-10)
  if (compareEmission) {
    check(sprintf("%s: max |emission params, EM - thmm1|", label),
          max(abs(packageParams(th$stPar, type) - fitted$emissionParameters)), emissionLimit)
  }
  invisible(list(th = th, em = em))
}

tt <- seq(0, 1, length.out = nPoints)
sinusoids <- t(sapply(0:(K - 1), function(j) sin(2 * pi * (tt + 0.2 * j))))
A1 <- matrix(0.09, K, K) + 0.55 * diag(K)     # the paper's A1: stay 0.64, move 0.09


# ===========================================================================
# 1. Structure: geometric-dwell HSMM == HMM
# ===========================================================================

cat("\n1. Geometric-dwell HSMM is the HMM with transition matrix A\n")

set.seed(2)
Arandom <- thmm1$rdirichlet(K, rep(2, K))
etaRandom <- as.vector(thmm1$rdirichlet(1, rep(2, K)))
geoModel <- thmmToHsmm(Arandom, etaRandom, cbind(seq(-4, 4, length.out = K), 1), "normal")
geoHmm <- hsmmToHmm(geoModel, initialType = "aggregateStart")

check("aggregate sizes - 1 (every geo state is one substate)", max(abs(geoHmm$aggregateSizes - 1)), 0)
check("max |hsmmToHmm()$gamma - A|", max(abs(geoHmm$gamma - Arandom)), 1e-14)
check("max |hsmmToHmm()$delta - eta|", max(abs(geoHmm$delta - etaRandom)), 1e-14)
check("max |A -> (omega, p) -> A|", max(abs(hsmmToThmmA(geoModel) - Arandom)), 1e-14)
check("max |dwell pmf - geo(1 - A_ii)| (d = 1..30)",
      max(vapply(seq_len(K), function(i) {
        max(abs(getModelDwellPmf(geoModel, i, 30) - dgeom(0:29, 1 - Arandom[i, i])))
      }, numeric(1))), 1e-12)

# Simulated chains: empirical transition matrices of simulateHsmm() and of
# thmm1.r's rMarkov() estimate the same A. With 1e5 steps the binomial sd of
# each entry is at most 0.5 / sqrt(1e5 / K) = 0.0035, so 0.02 is > 5 sd.
nChain <- 1e5
empiricalA <- function(states) {
  counts <- table(factor(head(states, -1), 1:K), factor(tail(states, -1), 1:K))
  unclass(counts / rowSums(counts))
}
pathPackage <- simulateHsmm(geoModel, nChain, initialType = "aggregateStart", seed = 3)$states
set.seed(4)
pathThmm <- thmm1$rMarkov(nChain, etaRandom, Arandom)
check("simulateHsmm(): max |empirical A - A|", max(abs(empiricalA(pathPackage) - Arandom)), 0.02)
check("thmm1 rMarkov(): max |empirical A - A|", max(abs(empiricalA(pathThmm) - Arandom)), 0.02)
check("max |empirical A, package - thmm1|", max(abs(empiricalA(pathPackage) - empiricalA(pathThmm))), 0.03)


# ===========================================================================
# 2. Fixed parameters: likelihood, posteriors, Viterbi
# ===========================================================================

cat("\n2. Likelihood, posteriors and Viterbi at fixed parameters\n")

compareFixed <- function(label, data, stPar, A, eta, type, norm = "W21", packageEmission = FALSE,
                         control = NULL) {
  th <- thmm1$thmm.init(data, stPar, A, eta, type = type)
  th$prInit <- log(th$prInit)
  th$prTran <- log(th$prTran)
  th <- setEmission(th, type, norm, packageEmission, control)
  th <- thmm1$thmm.alfa(th, TRUE, type)
  th <- thmm1$thmm.beta(th, TRUE, type)
  th <- thmm1$thmm.viterbi(th, type)
  logLikThmm <- thmm1$sumLogs.vec(th$alfa[nrow(th$alfa), ])
  posteriorThmm <- exp(th$alfa + th$beta - logLikThmm)

  model <- thmmToHsmm(A, eta, packageParams(stPar, type), packageType(type), control)
  x <- packageData(data, type)
  eStep <- hsmmEStep(model, x, stationary = FALSE)
  viterbi <- hsmmViterbiEM(model, x, stationary = FALSE)

  check(paste(label, ": |log L, package - thmm1|"), abs(eStep$logLik - logLikThmm), 1e-8)
  check(paste(label, ": max |state posterior, package - thmm1|"),
        max(abs(eStep$stateWeights - posteriorThmm)), 1e-8)
  check(paste(label, ": Viterbi paths differ at # times"), sum(viterbi$states != th$bestSeq), 0)
  check(paste(label, ": |Viterbi log-prob, package - thmm1|"), abs(viterbi$logProb - th$bestProb), 1e-8)
}

set.seed(137)
simNonpar <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, type = "nonpar",
                            len = nPoints, drift = sinusoids)
startMeans <- simNonpar$data[sample(nCurves, K), ]
A0 <- thmm1$rdirichlet(K, rep(6, K))
eta0 <- as.vector(thmm1$rdirichlet(1, rep(6, K)))

# At the truth (well separated: long sharply-peaked Viterbi path) and at a
# poor random start (diffuse posteriors).
for (norm in c("L2", "W21", "W22")) {
  compareFixed(sprintf("nonpar %s, truth", norm), simNonpar$data, sinusoids, A1, rep(1 / K, K),
               "nonpar", norm = norm, control = list(norm = norm))
  compareFixed(sprintf("nonpar %s, start", norm), simNonpar$data, startMeans, A0, eta0,
               "nonpar", norm = norm, control = list(norm = norm))
}

bmDrift <- c(-8, -4, 0, 4, 8)
set.seed(21)
simBm <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, drift = bmDrift,
                        type = "BMWD", len = nPoints)
compareFixed("BMWD, truth (package emission)", simBm$data, bmDrift, A1, rep(1 / K, K), "BMWD",
             packageEmission = TRUE)
compareFixed("BMWD, start (package emission)", simBm$data, c(-6, -2, 0, 2, 6), A0, eta0, "BMWD",
             packageEmission = TRUE)


# ===========================================================================
# 3. Baum-Welch iterates from a common start
# ===========================================================================

cat("\n3. Baum-Welch iterates, EM (geo dwell) vs thmm1.r, common start\n")

# 3a. Nonparametric mean curves: thmm1.r entirely unmodified.
for (norm in c("L2", "W21", "W22")) {
  compareIterates(sprintf("nonpar %s", norm), simNonpar$data, startMeans, A0, eta0, "nonpar",
                  norm = norm, control = list(norm = norm))
}

# 3b. Two curves per observation (thmm1.r's nonparMulti, W21 summed over curves).
multiTruth <- buildHsmm(
  {o <- A1 / (1 - diag(A1)); diag(o) <- 0; o}, lapply(1 - diag(A1), function(q) c(prob = q)),
  rep("geo", K), cbind(sinusoids, sinusoids[c(3:K, 1:2), ]), "nonpar", delta = rep(1 / K, K),
  stationary = FALSE, emissionControl = list(grid = tt, norm = "W21")
)
multiData <- simulateHsmm(multiTruth, nCurves, initialType = "aggregateStart", seed = 9)$x
set.seed(10)
multiStart <- multiData[sample(nCurves, K), , ]
compareIterates("nonparMulti (2 curves) W21", multiData, multiStart, A0, eta0, "nonparMulti",
                control = list(norm = "W21"))

# 3c. Brownian motion with drift, 1-d and 2-d: thmm1.r's forward-backward and
# re-estimation (endpoint-difference drift M-step), package emission matrix.
compareIterates("BMWD (package emission)", simBm$data, c(-6, -2, 0, 2, 6), A0, eta0, "BMWD",
                packageEmission = TRUE)

drift2d <- 3 * rbind(cos(2 * pi * (1:K) / K), sin(2 * pi * (1:K) / K))
set.seed(31)
simBm2d <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, drift = drift2d,
                          type = "BM2d", len = nPoints)
compareIterates("BM2d (package emission)", simBm2d$data, 0.7 * drift2d[, c(2:K, 1)], A0, eta0, "BM2d",
                packageEmission = TRUE)

# 3d. One iteration where only the dynamics M-step is exactly comparable.
set.seed(41)
ouPar <- rbind(c(4, 4, 8, 2, 20) * c(-2, 0, 4, 2, 1), c(4, 4, 8, 2, 20))
simOu <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, drift = ouPar[1, ] / ouPar[2, ],
                        theta = ouPar[2, ], type = "OU", len = nPoints)
compareIterates("OU, 1 iteration (package emission)", simOu$data, ouPar * 0.8, A0, eta0, "OU",
                iterations = 1, packageEmission = TRUE, compareEmission = FALSE)

# fBM (hurst 0.8): thmm1.r's drift M-step is coef * sum_k (k / L)^ups dO_k,
# the package's is int g dO / int g^2 with cell-averaged g = kappa tau^ups
# on the observation grid. Same continuum estimator, different quadrature,
# so the drifts agree only to O(1 / L).
set.seed(51)
fbmDrift <- c(-10, -6, -2, 0, 2)
simFbm <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, drift = fbmDrift,
                         type = "BMWD", len = nPoints, hurst = 0.8)
fbm <- compareIterates("fBM hurst 0.8, 1 iteration (package emission)", simFbm$data, fbmDrift + 1, A0, eta0,
                       "BMWD", iterations = 1, hurst = 0.8, packageEmission = TRUE,
                       control = list(hurst = 0.8), compareEmission = FALSE)
check("fBM hurst 0.8: max relative |drift, EM - thmm1| (quadrature)",
      max(abs(fbm$th$stPar - fbm$em$fittedModel$emissionParameters[, 1]) / pmax(abs(fbm$th$stPar), 1)), 0.05)


# ===========================================================================
# 4. Full fits to convergence on thmm1.r-simulated data
# ===========================================================================

cat("\n4. Fits to convergence on data simulated from a THMM (geometric truth)\n")

# 4a. thmm.BaumWelch() itself (its own stopping rule) vs fitHsmmEM, nonpar L2.
# (Not W21: thmm1.r's discrete W21 is not scaled by the grid spacing, so on
# 100 points the differenced Brownian-bridge noise swamps the signal and even
# the true parameters decode at ARI ~ 0; the iterates still match, see 3a.)
invisible(capture.output(
  bw <- thmm1$thmm.BaumWelch(simNonpar$data, K, state_mean = startMeans, prob.trans = A0, prob.init = eta0,
                             maxIter = 2000, minIter = 10, verbose = 1, tol = 1e-12, type = "nonpar",
                             norm = "L2")
))
emNonpar <- suppressWarnings(fitHsmmEM(
  thmmToHsmm(A0, eta0, startMeans, "nonpar", list(norm = "L2")), simNonpar$data,
  stationary = FALSE, tol = 1e-13, maxIter = 2000
))
vNonpar <- hsmmViterbiEM(emNonpar$fittedModel, simNonpar$data)
# thmm.BaumWelch() returns alfa from before its last re-estimation, so its
# log-likelihood is recomputed at the returned parameters.
bwFinal <- thmm1$thmm.alfa(thmm1$thmm.bj.all(bw, log.scale = TRUE, type = "nonpar", norm = "L2"), TRUE, "nonpar")
check("nonpar L2 converged: |log L, EM - thmm.BaumWelch|",
      abs(emNonpar$logLik - thmm1$sumLogs.vec(bwFinal$alfa[nrow(bwFinal$alfa), ])), 1e-6)
check("nonpar L2 converged: max |A, EM - thmm.BaumWelch|",
      max(abs(hsmmToThmmA(emNonpar$fittedModel) - exp(bw$prTran))), 1e-5)
check("nonpar L2 converged: max |mean curves, EM - thmm.BaumWelch|",
      max(abs(emNonpar$fittedModel$emissionParameters - bw$stPar)), 1e-5)
check("nonpar L2 converged: Viterbi paths differ at # times", sum(vNonpar$states != bw$bestSeq), 0)
check("nonpar L2: 1 - ARI(EM Viterbi, truth)", 1 - mclust::adjustedRandIndex(vNonpar$states, simNonpar$ss), 0.1)

# 4b. Monte Carlo, BMWD (paper's low/medium separation drifts), geometric A1
# truth. thmm1.r uses its own E-step functional (sum |dO - c / n|)^2, the
# package the Girsanov functional c S - c^2 / 2; both start from the same
# drifts and a flat A. Only statistical agreement is expected here.
nReplicates <- 6
mc <- data.frame()
flat <- matrix(1 / K, K, K)
for (replicate in seq_len(nReplicates)) {
  set.seed(1000 + replicate)
  sim <- thmm1$sim.thmm(n = nCurves, prInit = rep(1 / K, K), prTran = A1, drift = bmDrift,
                        type = "BMWD", len = nPoints)
  start <- sort(bmDrift + rnorm(K, sd = 1))

  invisible(capture.output(
    bwBm <- thmm1$thmm.BaumWelch(sim$data, K, state_mean = start, prob.trans = flat, prob.init = rep(1 / K, K),
                                 maxIter = 500, minIter = 10, verbose = 1, tol = 1e-10, type = "BMWD")
  ))
  emBm <- suppressWarnings(fitHsmmEM(thmmToHsmm(flat, rep(1 / K, K), matrix(start), "bmwd"), sim$data,
                                     stationary = FALSE, tol = 1e-10))
  vBm <- hsmmViterbiEM(emBm$fittedModel, sim$data)

  mc <- rbind(mc, data.frame(
    replicate = replicate,
    ariPackage = mclust::adjustedRandIndex(vBm$states, sim$ss),
    ariThmm1 = mclust::adjustedRandIndex(bwBm$bestSeq, sim$ss),
    ariBetween = mclust::adjustedRandIndex(vBm$states, bwBm$bestSeq),
    driftErrPackage = max(abs(emBm$fittedModel$emissionParameters[, 1] - bmDrift)),
    driftErrThmm1 = max(abs(bwBm$stPar - bmDrift)),
    driftBetween = max(abs(emBm$fittedModel$emissionParameters[, 1] - bwBm$stPar)),
    stayErrPackage = max(abs(diag(hsmmToThmmA(emBm$fittedModel)) - diag(A1))),
    stayErrThmm1 = max(abs(diag(exp(bwBm$prTran)) - diag(A1)))
  ))
}
cat("\n  BMWD Monte Carlo (each code with its own functional):\n")
print(format(mc, digits = 3), row.names = FALSE)
cat("\n")

check("BMWD MC: 1 - mean ARI(package, truth)", 1 - mean(mc$ariPackage), 0.15)
check("BMWD MC: 1 - mean ARI(thmm1, truth)", 1 - mean(mc$ariThmm1), 0.25)
check("BMWD MC: 1 - mean ARI(package, thmm1)", 1 - mean(mc$ariBetween), 0.25)
check("BMWD MC: mean ARI, thmm1 - package (package not worse)", mean(mc$ariThmm1 - mc$ariPackage), 0.02)
check("BMWD MC: mean max |drift - truth|, package", mean(mc$driftErrPackage), 1)
check("BMWD MC: mean max |drift, package - thmm1|", mean(mc$driftBetween), 1)
check("BMWD MC: mean max |A_ii - truth|, package", mean(mc$stayErrPackage), 0.15)


# ---------------------------------------------------------------------------

cat(sprintf("\n%d / %d checks passed\n", sum(checks$pass), nrow(checks)))
if (!all(checks$pass)) quit(status = 1)
