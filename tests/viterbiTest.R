# Does a THSMM with geometric dwell times converge to the same model as the
# THMM? One synthetic data set, no replicates.
#
# Run from the package root:   Rscript tests/viterbiTest.R
#
# With geometric dwell d_i ~ geo(p_i) and embedded jumps omega (zero
# diagonal), an HSMM is exactly the HMM with transition matrix
#   A = diag(1 - p) + diag(p) omega.
# So if both codes evaluate the same emission functional and start from the
# same point, EM on the THSMM and thmm1.R's Baum-Welch should produce the same
# iterates, the same fitted model and the same Viterbi path.
#
# Emissions are nonparametric mean curves with the L2 norm,
#   log b_j(O) = -1/2 sum_k (O(tau_k) - h_j(tau_k))^2,
# the one curve family where thmm1.R and the package use the identical
# functional. (For BM / OU / fBM, thmm1.R's E-step functional differs from the
# paper's, which the package follows, so the two fits would differ for that
# reason alone.)

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
codeFiles <- c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R",
               "likelihood.R")
for (file in codeFiles) source(file)

thmmPath <- file.path("..", "archivedCode", "thmm1.r")
if (!file.exists(thmmPath)) thmmPath <- file.path("tests", "thmm1.r")
thmm1 <- new.env()
suppressMessages(source(thmmPath, local = thmm1))   # also attaches mclust


# ---------------------------------------------------------------------------
# One synthetic data set
# ---------------------------------------------------------------------------

K <- 3
nCurves <- 300
grid <- seq(0, 1, length.out = 100)

# Three well-separated mean curves.
meanCurves <- rbind(
  sin(2 * pi * grid),
  2 * grid - 1,
  cos(2 * pi * grid) + 1
)

# Geometric dwell times: leave probabilities 0.3, 0.2, 0.1 (mean dwells 3.3,
# 5, 10), and uniform jumps between the other two states.
leave <- c(0.3, 0.2, 0.1)
omega <- matrix(0.5, K, K)
diag(omega) <- 0

truth <- buildHsmm(omega, lapply(leave, function(p) c(prob = p)), rep("geo", K), meanCurves, "nonpar",
                   delta = rep(1 / K, K), stationary = FALSE, emissionControl = list(grid = grid, norm = "L2"))
sim <- simulateHsmm(truth, nCurves, initialType = "aggregateStart", seed = 1)
x <- sim$x

matplot(grid, t(x), type = "l", lty = 1, col = c("#0072B2", "#D55E00", "#009E73")[sim$states],
        xlab = expression(tau), ylab = expression(O(tau)), main = "Simulated curves, coloured by true state")


# ---------------------------------------------------------------------------
# Common starting point
# ---------------------------------------------------------------------------

start <- initializeHsmmEM(x, K, "geo", "nonpar", list(grid = grid, norm = "L2"), stationary = FALSE, seed = 1)$model

# The THMM transition matrix of a geometric-dwell model.
thmmMatrix <- function(model) {
  hazard <- vapply(model$sojournParameters, function(s) s[["prob"]], numeric(1))
  diag(1 - hazard) + hazard * model$tpm
}


# ---------------------------------------------------------------------------
# Fit the THMM (thmm1.R) and the THSMM with geometric dwell (package EM)
# ---------------------------------------------------------------------------

invisible(capture.output(
  thmm <- thmm1$thmm.BaumWelch(x, K, state_mean = start$emissionParameters, prob.trans = thmmMatrix(start),
                               prob.init = start$delta, maxIter = 5000, minIter = 10, tol = 1e-14,
                               verbose = 1, type = "nonpar", norm = "L2")
))
# thmm.BaumWelch() keeps alfa from before its last re-estimation, so the
# log-likelihood is recomputed at the parameters it returns.
thmmFinal <- thmm1$thmm.alfa(thmm1$thmm.bj.all(thmm, log.scale = TRUE, type = "nonpar", norm = "L2"), TRUE, "nonpar")
thmmLogLik <- thmm1$sumLogs.vec(thmmFinal$alfa[nrow(thmmFinal$alfa), ])

thsmm <- fitHsmmEM(start, x, stationary = FALSE, tol = 1e-14, maxIter = 5000)
thsmmModel <- thsmm$fittedModel
thsmmViterbi <- hsmmViterbiEM(thsmmModel, x)


# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------

ari <- c(
  `THMM vs truth` = adjustedRandIndex(thmm$bestSeq, sim$states),
  `THSMM vs truth` = adjustedRandIndex(thsmmViterbi$states, sim$states),
  `THMM vs THSMM` = adjustedRandIndex(thmm$bestSeq, thsmmViterbi$states)
)

differences <- c(
  `|log-likelihood|` = abs(thmmLogLik - thsmm$logLik),
  `max |A|` = max(abs(exp(thmm$prTran) - thmmMatrix(thsmmModel))),
  `max |initial distribution|` = max(abs(exp(thmm$prInit) - thsmmModel$delta)),
  `max |mean curves|` = max(abs(thmm$stPar - thsmmModel$emissionParameters)),
  `Viterbi paths differ at # curves` = sum(thmm$bestSeq != thsmmViterbi$states),
  `|Viterbi log-probability|` = abs(thmm$bestProb - thsmmViterbi$logProb)
)

cat("\nLog-likelihood: THMM", format(thmmLogLik, digits = 12), " THSMM", format(thsmm$logLik, digits = 12),
    sprintf(" (THSMM: %d EM iterations, converged = %s)\n", thsmm$iterations, thsmm$converged))

cat("\nFitted transition matrix A\n  THMM (thmm1.R):\n")
print(round(exp(thmm$prTran), 4))
cat("  THSMM, geometric dwell (A = diag(1 - p) + diag(p) omega):\n")
print(round(thmmMatrix(thsmmModel), 4))
cat("  truth:\n")
print(round(thmmMatrix(truth), 4))

cat("\nViterbi ARI\n")
print(round(ari, 6))

cat("\nTHMM - THSMM differences\n")
print(signif(differences, 3))

cat("\nTrue states (rows) vs THSMM Viterbi states (columns)\n")
print(table(truth = sim$states, thsmm = thsmmViterbi$states))

same <- differences[["Viterbi paths differ at # curves"]] == 0 && max(differences[-5]) < 1e-6
cat(sprintf("\n%s\n", if (same) {
  "The THSMM with geometric dwell times converged to the same model as the THMM."
} else {
  "The THSMM and the THMM did NOT converge to the same model (see the differences above)."
}))
