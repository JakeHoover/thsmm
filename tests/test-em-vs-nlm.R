# EM vs numerical MLE (nlm) for HSMMs with Normal emissions.
#
# Run from the package root:  Rscript tests/test-em-vs-nlm.R
#
# What is checked
#   1. E-step machinery: the sparse forward pass reproduces the dense
#      hsmmNegativeLogLik(), the closed-form stationary distribution matches
#      getStationaryDistribution(), and sum_r d_i(r) = sum_j N(i, j).
#   2. stationary = FALSE (exact EM): from the same jittered start, EM and nlm
#      reach the same log-likelihood and parameters, EM is monotone, and the
#      numerical gradient of the nlm objective at the EM solution is ~0.
#      Three models cover every dwell family: shiftp + np; nb + weibull + npg
#      (non-trivial omega); geo + pois + shiftnb.
#   3. stationary = TRUE (plug-in delta): EM lands close to, not exactly at,
#      the stationary MLE; polish = TRUE closes the gap.
#   4. Shift profiling: fitHsmmShiftGrid() picks the same shift with either
#      fitter.

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")

for (f in list.files(".", pattern = "\\.R$", full.names = TRUE)) {
  if (basename(f) != "test-em-vs-nlm.R") source(f)
}

results <- data.frame(test = character(), value = numeric(), limit = numeric(), pass = logical())

check <- function(test, value, limit) {
  pass <- is.finite(value) && value <= limit
  results[nrow(results) + 1L, ] <<- list(test, value, limit, pass)
  cat(sprintf("  [%s] %-58s %10.3e  (limit %.0e)\n", if (pass) "PASS" else "FAIL", test, value, limit))
}

numericGradient <- function(fn, par, step = 1e-5) {
  vapply(seq_along(par), function(k) {
    up <- par; down <- par
    up[k] <- up[k] + step; down[k] <- down[k] - step
    (fn(up) - fn(down)) / (2 * step)
  }, numeric(1))
}

# Distance between two fitted models on the constrained scale. Dwell
# distributions are compared as pmfs (via the expanded-HMM representation),
# which is what the likelihood actually sees; raw nb/weibull parameters can be
# weakly identified even when the pmfs agree.
modelDistance <- function(a, b, maxDwell = 80) {
  dwell <- max(vapply(seq_len(a$nStates), function(i) {
    sum(abs(getModelDwellPmf(a, i, maxDwell) - getModelDwellPmf(b, i, maxDwell)))
  }, numeric(1)))

  c(
    emission = max(abs(a$emissionParameters - b$emissionParameters)),
    tpm = max(abs(a$tpm - b$tpm)),
    delta = if (is.null(a$delta) || is.null(b$delta)) 0 else max(abs(a$delta - b$delta)),
    dwellL1 = dwell
  )
}


# ---------------------------------------------------------------------------
# Test models
# ---------------------------------------------------------------------------

modelA <- exampleHsmmModel()   # shiftp (shift fixed at 1) + np, K = 2

modelB <- buildHsmm(
  tpm = matrix(c(0, .7, .3,
                 .5, 0, .5,
                 .2, .8, 0), 3, byrow = TRUE),
  sojournParameters = list(
    c(size = 3, prob = .4, dwellSize = 20),
    c(shape = 1.8, scale = 5, dwellSize = 15),
    c(.2, .3, .2, .1, .25)
  ),
  sojournType = c("nb", "weibull", "npg"),
  emissionParameters = matrix(c(-2, 1, 1, .8, 4, 1.2), 3, byrow = TRUE),
  emissionType = "normal",
  delta = c(.5, .3, .2)
)

modelC <- buildHsmm(
  tpm = matrix(c(0, .6, .4,
                 .3, 0, .7,
                 .5, .5, 0), 3, byrow = TRUE),
  sojournParameters = list(
    c(prob = .3),
    c(lambda = 5, dwellSize = 20),
    c(size = 2, prob = .5, shift = 3, dwellSize = 20)
  ),
  sojournType = c("geo", "pois", "shiftnb"),
  emissionParameters = matrix(c(0, 1, 3, 1, 6, 1.5), 3, byrow = TRUE),
  emissionType = "normal",
  delta = c(.4, .4, .2)
)

models <- list(A = modelA, B = modelB, C = modelC)
nObs <- c(A = 1000, B = 1500, C = 1500)


# ---------------------------------------------------------------------------
# 1. E-step machinery
# ---------------------------------------------------------------------------

cat("\n1. E-step machinery\n")

for (name in names(models)) {
  model <- models[[name]]
  x <- simulateHsmm(model, 400, seed = 3)$x
  x[c(10, 50:52)] <- NA                                 # exercise missing data

  for (stationary in c(TRUE, FALSE)) {
    eStep <- hsmmEStep(model, x, stationary)
    reference <- -hsmmNegativeLogLik(
      hsmmConstrainedToWorking(model, stationary), x, model, stationary = stationary
    )
    deltaDense <- hsmmToHmm(model, if (stationary) "stationaryHmm" else "aggregateStart")$delta
    tag <- sprintf("%s, stationary=%s", name, stationary)

    check(paste(tag, ": |logLik sparse - dense|"), abs(eStep$logLik - reference), 1e-8)
    check(paste(tag, ": |delta closed form - dense|"),
          max(abs(getExpandedInitial(model, eStep$chain, stationary) - deltaDense)), 1e-10)
    check(paste(tag, ": |sum_r d - sum_j N|"),
          max(abs(rowsum(eStep$departures, eStep$chain$aggregate)[, 1] - rowSums(eStep$transitions))),
          1e-8)
  }
}


# ---------------------------------------------------------------------------
# 2. Exact EM (stationary = FALSE) vs nlm
# ---------------------------------------------------------------------------

cat("\n2. EM vs nlm, stationary = FALSE\n")

timing <- data.frame()

for (name in names(models)) {
  trueModel <- models[[name]]
  x <- simulateHsmm(trueModel, nObs[[name]], initialType = "aggregateStart", seed = 21)$x
  set.seed(5)
  start <- jitterHsmmStart(trueModel, x, emissionScale = 0.3, dwellScale = 0.3)

  emTime <- system.time(em <- fitHsmmEM(start, x, stationary = FALSE))[["elapsed"]]
  nlmTime <- system.time(
    nl <- suppressWarnings(fitHsmm(start, x, stationary = FALSE, iterlim = 1000))
  )[["elapsed"]]

  distance <- modelDistance(em$fittedModel, nl$fittedModel)

  objective <- function(par) hsmmNegativeLogLik(par, x, start, stationary = FALSE)
  gradientAtEm <- numericGradient(objective, hsmmConstrainedToWorking(em$fittedModel, FALSE))

  cat(sprintf("  model %s: logLik EM %.6f | nlm %.6f\n", name, em$logLik, -nl$minimum))
  check(paste(name, ": EM - nlm logLik (EM should not be worse)"), -nl$minimum - em$logLik, 1e-4)
  check(paste(name, ": |EM - nlm| logLik"), abs(em$logLik + nl$minimum), 1e-3)
  check(paste(name, ": EM monotone (largest decrease)"), max(0, -diff(em$logLikTrace)), 1e-8)
  check(paste(name, ": max |grad nll| at EM solution"), max(abs(gradientAtEm)), 1e-2)
  check(paste(name, ": max |emission diff|"), distance[["emission"]], 1e-3)
  check(paste(name, ": max |omega diff|"), distance[["tpm"]], 1e-3)
  check(paste(name, ": max |delta diff|"), distance[["delta"]], 1e-2)
  check(paste(name, ": max L1 dwell-pmf diff"), distance[["dwellL1"]], 1e-2)

  timing <- rbind(timing, data.frame(
    model = name, emIterations = em$iterations, emSeconds = emTime,
    nlmIterations = nl$iterations, nlmSeconds = nlmTime, nlmCode = nl$code
  ))
}

cat("\n  Cost comparison:\n")
print(timing, row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. Stationary case: plug-in EM, then polish
# ---------------------------------------------------------------------------

cat("\n3. stationary = TRUE (plug-in delta)\n")

x <- simulateHsmm(modelB, 1500, seed = 21)$x
set.seed(5)
start <- jitterHsmmStart(modelB, x, emissionScale = 0.3, dwellScale = 0.3)

emStat <- fitHsmmEM(start, x, stationary = TRUE)
nlStat <- suppressWarnings(fitHsmm(start, x, stationary = TRUE, iterlim = 1000))
polished <- suppressWarnings(fitHsmmEM(start, x, stationary = TRUE, polish = TRUE))

cat(sprintf("  logLik  plug-in EM %.6f | nlm %.6f | EM + polish %.6f\n",
            emStat$logLik, -nlStat$minimum, -polished$minimum))
check("plug-in EM gap to stationary MLE (small, not ~0)", -nlStat$minimum - emStat$logLik, 0.5)
check("|EM + polish - nlm| logLik", abs(polished$minimum - nlStat$minimum), 1e-3)
check("polish nlm iterations (should be few)", polished$iterations, 100)


# ---------------------------------------------------------------------------
# 4. Shift profiling with either fitter
# ---------------------------------------------------------------------------

cat("\n4. Shift grid (model C, state 3 true shift = 3)\n")

x <- simulateHsmm(modelC, 1500, initialType = "aggregateStart", seed = 21)$x
gridArgs <- list(model = modelC, dat = x, shiftValues = list(NULL, NULL, 1:5), stationary = FALSE)

gridEm <- do.call(fitHsmmShiftGrid, c(gridArgs, fitFunction = fitHsmmEM))
gridNlm <- suppressWarnings(do.call(fitHsmmShiftGrid, gridArgs))

print(data.frame(shift = 1:5, nllEM = gridEm$shiftGrid$negativeLogLik,
                 nllNlm = gridNlm$shiftGrid$negativeLogLik))
# At a wrong shift the best NB can sit on the edge of the family (size -> Inf,
# i.e. the Poisson limit): the supremum is not attained and both fitters drift
# along a flat ridge, nlm usually stopping first with code 3. So away from the
# best shift only "EM is not worse" is meaningful; at the best shift (interior
# MLE) the two must agree.
bestRow <- which.min(gridNlm$shiftGrid$negativeLogLik)
check("same best shift", abs(gridEm$bestShift[3] - gridNlm$bestShift[3]), 0)
check("EM nll - nlm nll, worst grid point (EM not worse)",
      max(gridEm$shiftGrid$negativeLogLik - gridNlm$shiftGrid$negativeLogLik), 1e-4)
check("|nll EM - nll nlm| at the best shift",
      abs(gridEm$shiftGrid$negativeLogLik[bestRow] - gridNlm$shiftGrid$negativeLogLik[bestRow]), 1e-3)


# ---------------------------------------------------------------------------

cat(sprintf("\n%d / %d checks passed\n", sum(results$pass), nrow(results)))
if (!all(results$pass)) quit(status = 1)
