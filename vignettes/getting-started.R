# Getting started: simulate, initialize, and fit a hidden semi-Markov model.
#
# Run from the package root:   Rscript vignettes/getting-started.R
# (about 30 seconds; plots go to Rplots.pdf), or better, step through it line by
# line in RStudio and look at each result.
#
# A hidden semi-Markov model (HSMM) has three ingredients:
#   * an embedded transition matrix `tpm` (omega): where the chain goes when it
#     leaves a state. Its diagonal is zero -- staying is handled by the dwell time.
#   * a dwell-time (sojourn) distribution per state: how long the chain stays.
#   * an emission distribution per state: what is observed while in a state.
#
# Part 1 does this for ordinary numbers (normal emissions). Part 2 does the
# same for curves -- the topological HMM (THMM) -- where each observation is a
# whole function, here Brownian motion with a state-dependent drift, and dwell
# times are geometric. Part 3 gives the curve states real dwell-time
# distributions -- the topological HSMM (THSMM) -- and compares the fits.

setwd("~/Documents/statistics/functionalHiddenModels/thsmm/code/V1")
for (file in c("dwell.R", "em.R", "initialize.R", "simulate.R", "transform.R",
               "emissionFamilies.R", "fit.R", "model.R", "viterbi.R", "likelihood.R")) {
  source(file)
}


# ===========================================================================
# Part 1: a 3-state HSMM with normal emissions
# ===========================================================================

# ---------------------------------------------------------------------------
# 1.1 Write down a "true" model
# ---------------------------------------------------------------------------

# Embedded transition matrix: rows sum to one, zero diagonal.
tpm <- matrix(c(0.0, 0.7, 0.3,
                0.5, 0.0, 0.5,
                0.4, 0.6, 0.0), nrow = 3, byrow = TRUE)

# One dwell distribution per state. Each family has its own parameter vector;
# `dwellSize` is how many dwell times are modelled exactly before a geometric
# tail takes over (make it comfortably larger than the typical dwell).
#   "pois"  c(lambda, dwellSize)          dwell = 1 + Poisson(lambda)
#   "nb"    c(size, prob, dwellSize)      dwell = 1 + negative binomial
#   "geo"   c(prob)                       geometric: this is an ordinary HMM
sojournType <- c("pois", "nb", "geo")
sojournParameters <- list(
  c(lambda = 6, dwellSize = 25),
  c(size = 3, prob = 0.3, dwellSize = 30),
  c(prob = 0.25)
)

# Emissions: one row per state, columns (mean, sd) for "normal".
emissionParameters <- cbind(mean = c(-3, 0, 4), sd = c(1, 0.8, 1.2))

trueModel <- buildHsmm(
  tpm = tpm,
  sojournParameters = sojournParameters,
  sojournType = sojournType,
  emissionParameters = emissionParameters,
  emissionType = "normal",
  delta = c(1, 1, 1) / 3,     # distribution of the first state
  stationary = FALSE          # start from `delta` rather than the stationary law
)


# ---------------------------------------------------------------------------
# 1.2 Simulate data
# ---------------------------------------------------------------------------

sim <- simulateHsmm(trueModel, nObs = 2000, initialType = "aggregateStart", seed = 1)

x <- sim$x                    # the observations
head(x)
table(sim$states)             # the hidden states (unknown in practice)
head(getStateRuns(sim$states))  # realized dwell times, one row per visit

plot(x[1:300], type = "l", col = "grey50", xlab = "time", ylab = "x",
     main = "First 300 observations, coloured by hidden state")
points(x[1:300], col = sim$states[1:300] + 1, pch = 19, cex = 0.5)


# ---------------------------------------------------------------------------
# 1.3 Initial model parameters
# ---------------------------------------------------------------------------

# In practice we do not know the truth. initializeHsmmEM() builds a starting
# model from the data alone: k-means on the observations gives emission
# parameters and a rough state labelling, from which it guesses the embedded
# transitions and mean dwell times. We only choose the number of states and
# the dwell family for each state.
init <- initializeHsmmEM(
  x,
  nStates = 3,
  sojournType = c("pois", "nb", "geo"),
  emissionType = "normal",
  dwellSize = 30,
  stationary = FALSE,
  seed = 1
)

startModel <- init$model
startModel$emissionParameters   # starting means and sds (states ordered by mean)
startModel$tpm
startModel$sojournParameters

# You can also build a starting model by hand with buildHsmm(), exactly as in
# 1.1 -- useful when you have prior knowledge about the states.


# ---------------------------------------------------------------------------
# 1.4 Fit the model
# ---------------------------------------------------------------------------

# fitHsmmEM() runs the EM algorithm. (fitHsmm() maximizes the same likelihood
# numerically with nlm(); EM is usually faster and is the only option for
# curve data.)
fit <- fitHsmmEM(startModel, x, stationary = FALSE)

fit$converged
fit$iterations
fit$logLik

fitted <- fit$fittedModel
round(fitted$emissionParameters, 2)   # compare with emissionParameters above
round(fitted$tpm, 2)                  # compare with tpm above

# EM never decreases the log-likelihood:
plot(fit$logLikTrace, type = "l", xlab = "EM iteration", ylab = "log-likelihood")


# ---------------------------------------------------------------------------
# 1.5 Look at the results
# ---------------------------------------------------------------------------

# Dwell-time distributions: fitted (points) vs truth (bars).
par(mfrow = c(1, 3))
for (i in 1:3) {
  plot(1:25, getModelDwellPmf(trueModel, i, 25), type = "h", lwd = 4, col = "grey70",
       xlab = "dwell time", ylab = "probability", main = paste("State", i, "-", sojournType[i]))
  points(1:25, getModelDwellPmf(fitted, i, 25), pch = 19, col = "firebrick")
}
par(mfrow = c(1, 1))

# Most likely hidden state sequence (Viterbi) ...
decoded <- hsmmViterbiEM(fitted, x)
table(decoded = decoded$states, truth = sim$states)
mean(decoded$states == sim$states)    # share of time points decoded correctly

# ... and the probability of each state at each time.
posterior <- hsmmStatePosterior(fitted, x)
head(round(posterior, 3))


# ===========================================================================
# Part 2: curves (topological HMM), Brownian motion with drift
# ===========================================================================

# Each observation is now a curve on a grid over [0, 1]. With emissionType
# "bmwd", state j's curves are Brownian motion with drift c_j. Other curve
# families: "ou" (Ornstein-Uhlenbeck) and "nonpar" (state mean curves).
# Settings for the family go in `emissionControl`.
#
# Geometric dwell times ("geo", c(prob)) make the HSMM an ordinary HMM: each
# step the chain leaves state i with probability prob_i, whatever the time
# already spent there. This is the THMM of Kashlak, Loliencar and Heo (2023).

grid <- seq(0, 1, length.out = 100)
K <- 3

curveTruth <- buildHsmm(
  tpm = matrix(c(0, 0.5, 0.5,
                 0.5, 0, 0.5,
                 0.5, 0.5, 0), 3, byrow = TRUE),
  sojournParameters = list(c(prob = 0.25), c(prob = 0.15), c(prob = 0.35)),
  sojournType = rep("geo", K),     # mean dwell 1 / prob = 4, 6.7, 2.9 curves
  emissionParameters = matrix(c(-4, 0, 4), ncol = 1),   # one drift per state
  emissionType = "bmwd",
  emissionControl = list(grid = grid),
  delta = c(1, 0, 0),
  stationary = FALSE
)

curveSim <- simulateHsmm(curveTruth, nObs = 300, initialType = "aggregateStart", seed = 2)
curves <- curveSim$x          # a 300 x 100 matrix: one curve per row
dim(curves)

matplot(grid, t(curves[1:60, ]), type = "l", lty = 1, col = curveSim$states[1:60] + 1,
        xlab = "tau", ylab = "O(tau)", main = "First 60 curves, coloured by hidden state")

# Initialize and fit exactly as before; only the emission type and control change.
curveInit <- initializeHsmmEM(curves, nStates = K, sojournType = "geo", emissionType = "bmwd",
                              emissionControl = list(grid = grid),
                              stationary = FALSE, seed = 2)
curveFit <- fitHsmmEM(curveInit$model, curves, stationary = FALSE)

round(curveFit$fittedModel$emissionParameters, 2)   # truth: -4, 0, 4
round(unlist(curveFit$fittedModel$sojournParameters), 2)   # truth: 0.25, 0.15, 0.35
curveDecoded <- hsmmViterbiEM(curveFit$fittedModel, curves)
mean(curveDecoded$states == curveSim$states)


# ===========================================================================
# Part 3: topological HSMM -- curves with dwell-time distributions
# ===========================================================================

# Geometric dwell times are memoryless and their most likely value is always
# 1. Real regimes often last a minimum time and then end at a typical length.
# Here the truth uses shifted-Poisson dwell times, c(lambda, shift, dwellSize):
# dwell = shift + Poisson(lambda), so state 1 always lasts at least 4 curves.
# The drifts are close together (-1.5, 0, 1.5), so single curves are hard
# to classify and knowing how long states last really helps.

dwellTruth <- buildHsmm(
  tpm = matrix(c(0, 0.5, 0.5,
                 0.5, 0, 0.5,
                 0.5, 0.5, 0), 3, byrow = TRUE),
  sojournParameters = list(c(lambda = 2, shift = 4, dwellSize = 20),   # mean dwell 6
                           c(lambda = 1, shift = 3, dwellSize = 20),   # mean dwell 4
                           c(lambda = 3, shift = 2, dwellSize = 20)),  # mean dwell 5
  sojournType = rep("shiftp", K),
  emissionParameters = matrix(c(-1.5, 0, 1.5), ncol = 1),
  emissionType = "bmwd",
  emissionControl = list(grid = grid),
  delta = c(1, 0, 0),
  stationary = FALSE
)

dwellSim <- simulateHsmm(dwellTruth, nObs = 500, initialType = "aggregateStart", seed = 3)
dwellCurves <- dwellSim$x
table(getStateRuns(dwellSim$states)$dwellTime)   # no visit is shorter than 2 curves


# ---------------------------------------------------------------------------
# 3.1 Fit three dwell-time models to the same curves
# ---------------------------------------------------------------------------

# Only `sojournType` changes between fits; the emission start is the same
# because initializeHsmmEM() is deterministic given the seed. tol = 1e-7
# keeps the tutorial quick (the default, 1e-10, is stricter).

# (a) Geometric dwell: the THMM, which ignores how long states last.
geoStart <- initializeHsmmEM(dwellCurves, nStates = K, sojournType = "geo", emissionType = "bmwd",
                             emissionControl = list(grid = grid), stationary = FALSE, seed = 3)
geoFit <- fitHsmmEM(geoStart$model, dwellCurves, stationary = FALSE, tol = 1e-7)

# (b) Nonparametric dwell, "npg": free probabilities for dwell times 1 to 9
# plus a geometric tail beyond them. Assumes no shape, but has more parameters.
npgStart <- initializeHsmmEM(dwellCurves, nStates = K, sojournType = "npg", emissionType = "bmwd",
                             emissionControl = list(grid = grid), dwellSize = 10,
                             stationary = FALSE, seed = 3)
npgFit <- fitHsmmEM(npgStart$model, dwellCurves, stationary = FALSE, tol = 1e-7)

# (c) Shifted Poisson, the true family. EM estimates lambda but keeps each
# shift fixed (it is an integer), so fitHsmmShiftGrid() refits the model for
# every combination of candidate shifts and keeps the best one. Here that is
# shifts 2, 3, 4 for each state: 3^3 = 27 EM fits (about 20 seconds).
shiftStart <- initializeHsmmEM(dwellCurves, nStates = K, sojournType = "shiftp", emissionType = "bmwd",
                               emissionControl = list(grid = grid), dwellSize = 20,
                               stationary = FALSE, seed = 3)
shiftFit <- fitHsmmShiftGrid(shiftStart$model, dwellCurves, shiftValues = list(2:4, 2:4, 2:4),
                             stationary = FALSE, fitFunction = fitHsmmEM, tol = 1e-7)
shiftFit$bestShift                 # truth: 4, 3, 2
head(shiftFit$shiftGrid[order(shiftFit$shiftGrid$negativeLogLik), ])


# ---------------------------------------------------------------------------
# 3.2 Compare the fits
# ---------------------------------------------------------------------------

# Number of free parameters, for BIC: emissions, dwell parameters (shifts are
# not counted), embedded transitions (K - 2 free entries per row) and delta.
countParameters <- function(model) {
  dwell <- sum(mapply(getSojournParamCount, model$sojournParameters, model$sojournType))
  length(model$emissionParameters) + dwell + K * (K - 2) + (K - 1)
}

dwellFits <- list(`THMM (geo)` = geoFit, `THSMM (npg)` = npgFit, `THSMM (shiftp)` = shiftFit)

comparison <- do.call(rbind, lapply(names(dwellFits), function(name) {
  model <- dwellFits[[name]]$fittedModel
  decoded <- hsmmViterbiEM(model, dwellCurves)
  data.frame(
    model = name,
    logLik = dwellFits[[name]]$logLik,
    nPar = countParameters(model),
    BIC = -2 * dwellFits[[name]]$logLik + countParameters(model) * log(nrow(dwellCurves)),
    accuracy = mean(decoded$states == dwellSim$states),
    drifts = paste(round(model$emissionParameters[, 1], 2), collapse = ", ")
  )
}))
comparison

# What to look for:
#   * npg has the highest log-likelihood, but it also has the most
#     parameters. BIC (lower is better) picks the true shifted-Poisson model.
#   * Decoding accuracy rises from about 61% (geo) to 84% (shiftp). The THMM
#     also badly misplaces the middle drift (truth 0).
#   * The grid may swap shifts between states whose mean dwells are close;
#     the dwell *distributions* below are what matter for decoding.

# Fitted dwell distributions vs the truth. States are matched by drift: every
# initializer orders the states by drift, as the truth does.
par(mfrow = c(1, 3))
for (i in 1:K) {
  plot(1:15, getModelDwellPmf(dwellTruth, i, 15), type = "h", lwd = 6, col = "grey75",
       ylim = c(0, 0.7), xlab = "dwell time (curves)", ylab = "probability",
       main = sprintf("State %d (drift %.1f)", i, dwellTruth$emissionParameters[i, 1]))
  lines(1:15, getModelDwellPmf(geoFit$fittedModel, i, 15), type = "b", pch = 19, col = "steelblue")
  lines(1:15, getModelDwellPmf(npgFit$fittedModel, i, 15), type = "b", pch = 19, col = "darkorange")
  lines(1:15, getModelDwellPmf(shiftFit$fittedModel, i, 15), type = "b", pch = 19, col = "forestgreen")
  if (i == 1) legend("topright", c("truth", "geo", "npg", "shiftp"), bty = "n",
                     col = c("grey75", "steelblue", "darkorange", "forestgreen"), lwd = c(6, 1, 1, 1))
}
par(mfrow = c(1, 1))


# ---------------------------------------------------------------------------
# Where next
# ---------------------------------------------------------------------------
#   * Other dwell families: "pois", "nb", "shiftnb", "weibull", "np"
#     (see ?getDwellPmf in dwell.R), and they can be mixed across states.
#   * tests/testsMaster.R and tests/paperSims.R are the simulation studies
#     used for the paper results.
