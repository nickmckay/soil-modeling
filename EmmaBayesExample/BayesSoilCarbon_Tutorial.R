# =============================================================================
# Bayesian Estimation of Soil Carbon Model Parameters
# =============================================================================
#
# A teaching example: use MCMC to estimate decomposition rates (k) from
# radiocarbon (14C) and carbon stock observations from permafrost soils
# at the CiPEHR warming experiment in interior Alaska.
#
# We start with a simple 1-pool model (2 parameters: k and input), then
# move to a 3-pool series model (6 parameters). This lets you see how
# Bayesian inference works before tackling a harder problem.
#
# The data come from soil cores measured in 2009 and 2022, with three soil
# layers: surface organic, deep organic, and mineral. For the 1-pool model
# we aggregate these into a single "bulk soil" pool.
#
# Dependencies: SoilR, BayesianTools, deSolve, ggplot2, dplyr, tidyr, patchwork
#
# Install if needed:
#   install.packages(c("BayesianTools", "ggplot2", "dplyr", "tidyr", "patchwork"))
#   install.packages("SoilR", repos = "http://R-Forge.R-project.org")
# =============================================================================

library(SoilR)
library(BayesianTools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(42)  # for reproducibility of MCMC

# =============================================================================
# SECTION 0: Load data
# =============================================================================
# The .RData file contains pre-processed observations from the CiPEHR
# experiment, plus an atmospheric 14C record (the "bomb curve").
#
# Key objects loaded:
#   datComb  - soil 14C and C stock observations (2009 initial + 2022 treatments)
#   atmIn    - atmospheric 14C record from 1900-2022 (columns: YEAR, Atmosphere)
#              combines pre-bomb (IntCal) + post-bomb (Hua) + site measurements

# Locate the data file relative to this script.
# Works whether you source() this script or run it with Rscript.
.script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),                          # works with source()
  error = function(e) {
    if (exists("argv")) dirname(argv[1])                 # Rscript on some systems
    else if (length(commandArgs(trailingOnly = FALSE)) > 0) {
      arg <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("^--file=", arg, value = TRUE)
      if (length(file_arg)) dirname(sub("^--file=", "", file_arg))
      else "."
    } else "."
  }
)
load(file.path(.script_dir, "soilRmcmc.RData"))

# Take a look at what we're working with
cat("=== Soil observations ===\n")
print(as.data.frame(datComb))

# --- Atmospheric 14C: the "bomb spike" ---
# Nuclear weapons testing in the 1950s-60s nearly doubled atmospheric 14C.
# This pulse propagates into soils at a rate governed by decomposition —
# fast-cycling pools track the atmosphere closely, slow pools barely respond.
# That's what makes 14C so powerful for estimating turnover times.

ggplot(subset(atmIn, YEAR >= 1940), aes(x = YEAR, y = Atmosphere)) +
  geom_line(linewidth = 1, color = "steelblue") +
  labs(title = "Atmospheric 14C: The Bomb Spike",
       subtitle = "Nuclear testing doubled atmospheric 14C in the 1960s; it has declined since",
       x = "Year", y = expression(Delta^14 * "C (per mil)")) +
  theme_minimal(base_size = 14)


# #############################################################################
#
#  PART 1: ONE-POOL MODEL  (2 parameters: k, input)
#
# #############################################################################
#
# Simplest case: treat the entire soil column as a single well-mixed pool.
# We aggregate the three layers (surface organic, deep organic, mineral)
# into bulk totals for C stock and a stock-weighted mean for 14C.
#
# Parameters to estimate:
#   k     = decomposition rate (yr^-1). Turnover time = 1/k.
#   input = carbon input rate (gC/m2/yr)

# =============================================================================
# STEP 1.1: Prepare the 1-pool observations
# =============================================================================
# Aggregate layered data into a single bulk pool per time point.
# C stocks: sum across layers. 14C: stock-weighted mean across layers.

# Use the Control treatment (ambient conditions)
dat_1pool <- datComb %>%
  filter(year == 2009 | treatment == "Control") %>%
  mutate(C_stock_g    = C_stock * 1000,           # convert kg/m2 -> gC/m2
         C_stock_sd_g = C_stock_sd * 1000) %>%
  group_by(year) %>%
  summarise(
    C14        = weighted.mean(C14, C_stock_g),    # stock-weighted mean 14C
    C_14_sd    = sqrt(sum((C_stock_g * C_14_sd)^2)) / sum(C_stock_g),
    C_stock    = sum(C_stock_g),
    C_stock_sd = sqrt(sum(C_stock_sd_g^2)),        # propagate uncertainty
    .groups = "drop"
  )

cat("\n=== 1-Pool bulk observations ===\n")
print(as.data.frame(dat_1pool))

# Initial conditions from 2009
C0_1pool <- dat_1pool$C_stock[dat_1pool$year == 2009]
F0_1pool <- dat_1pool$C14[dat_1pool$year == 2009]
years_1pool <- seq(2009, 2022)

cat("\nInitial C stock:", round(C0_1pool), "gC/m2")
cat("\nInitial Delta14C:", round(F0_1pool, 1), "per mil\n")

# =============================================================================
# STEP 1.2: Define the likelihood function
# =============================================================================
# The likelihood measures how well proposed parameter values explain the data.
# For each candidate (k, input), we:
#   1. Run the forward model from 2009 to 2022
#   2. Compare predictions at 2022 to observations
#   3. Return the log-likelihood (higher = better fit)

obs_2022_1pool <- dat_1pool %>% filter(year == 2022)

likelihood_1pool <- function(pars) {
  k     <- pars[1]
  input <- pars[2]

  # Run the 1-pool model forward from 2009 to 2022
  model <- tryCatch(
    OnepModel14(
      t = years_1pool,
      k = k,
      C0 = C0_1pool,
      F0_Delta14C = F0_1pool,
      In = input,
      inputFc = atmIn,
      pass = TRUE
    ),
    error = function(e) return(NULL)
  )

  if (is.null(model)) return(-1e10)

  # Get predictions at 2022 (last time step)
  pred_C14    <- getF14(model)[length(years_1pool), ]
  pred_Cstock <- getC(model)[length(years_1pool), ]

  # Log-likelihood: Gaussian errors with observed measurement uncertainty
  ll_C14    <- dnorm(pred_C14 - obs_2022_1pool$C14, 0,
                     obs_2022_1pool$C_14_sd, log = TRUE)
  ll_Cstock <- dnorm(pred_Cstock - obs_2022_1pool$C_stock, 0,
                     obs_2022_1pool$C_stock_sd, log = TRUE)

  return(sum(ll_C14, ll_Cstock, na.rm = TRUE))
}

# Quick sanity check
cat("\nTest likelihood at k=0.01, input=500:",
    likelihood_1pool(c(0.01, 500)), "\n")

# =============================================================================
# STEP 1.3: Define the prior distribution
# =============================================================================
# The prior encodes what we know (or assume) BEFORE seeing the data.
# A uniform prior says "the parameter is equally likely anywhere in this range."
#
# For a bulk permafrost soil profile, reasonable ranges might be:
#   k: 0.001 to 0.1 yr^-1 (turnover time 10 to 1000 years)
#   input: 50 to 1000 gC/m2/yr
#
# *** EXERCISE: Try changing these bounds and see how it affects the posterior! ***
#   - What happens with a very wide prior for k (e.g., 0.0001 to 1)?
#   - What happens if you constrain input tightly (e.g., 200 to 300)?
#   - What if the prior excludes the true region entirely?

prior_lower_1pool <- c(k = 0.001, input = 50)
prior_upper_1pool <- c(k = 0.1,   input = 1000)

prior_1pool <- createUniformPrior(
  lower = prior_lower_1pool,
  upper = prior_upper_1pool
)

# =============================================================================
# STEP 1.4: Visualize the prior
# =============================================================================
# Always look at your prior before running MCMC!

prior_samples_1pool <- data.frame(
  k     = runif(10000, prior_lower_1pool["k"],     prior_upper_1pool["k"]),
  input = runif(10000, prior_lower_1pool["input"], prior_upper_1pool["input"])
)

prior_long_1pool <- prior_samples_1pool %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c("k", "input"),
                            labels = c("k (yr^-1)", "input (gC/m2/yr)")))

prior_plot_1pool <- ggplot(prior_long_1pool, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "skyblue", color = "white", alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +
  labs(title = "Prior Distributions (1-Pool Model)",
       subtitle = "Uniform priors --every value in the range is equally likely a priori",
       x = "Parameter value", y = "Density") +
  theme_minimal(base_size = 14)

print(prior_plot_1pool)

# =============================================================================
# STEP 1.5: Run the MCMC
# =============================================================================
# The MCMC sampler explores parameter space, spending more time in regions
# of high posterior probability (good fit to data AND consistent with prior).
#
# DREAMzs runs multiple interacting chains --this helps it explore
# complicated posterior shapes. For a 2-parameter problem, 10,000
# iterations is usually plenty.

setup_1pool <- createBayesianSetup(
  likelihood = likelihood_1pool,
  prior      = prior_1pool,
  names      = c("k", "input")
)

settings_1pool <- list(iterations = 10000, message = FALSE)

cat("\nRunning 1-pool MCMC (this should take ~10-30 seconds)...\n")
mcmc_1pool <- runMCMC(
  bayesianSetup = setup_1pool,
  sampler       = "DREAMzs",
  settings      = settings_1pool
)
cat("Done!\n")

# =============================================================================
# STEP 1.6: Visualize the MCMC chains (trace plots)
# =============================================================================
# Trace plots show the path each chain took through parameter space.
# Good chains should:
#   - Look like "hairy caterpillars" (well-mixed, no trends)
#   - Overlap with each other (converged to same region)
#   - NOT show long drifts or get stuck in one place

raw_chains_1 <- getSample(mcmc_1pool, start = 1, coda = TRUE)

chain_df_1 <- do.call(rbind, lapply(seq_along(raw_chains_1), function(i) {
  ch <- as.data.frame(as.matrix(raw_chains_1[[i]]))
  colnames(ch) <- c("k", "input")
  ch$iteration <- seq_len(nrow(ch))
  ch$chain <- factor(i)
  ch
}))

chain_long_1 <- chain_df_1 %>%
  pivot_longer(cols = c("k", "input"),
               names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c("k", "input"),
                            labels = c("k (yr^-1)", "input (gC/m2/yr)")))

trace_1pool <- ggplot(chain_long_1, aes(x = iteration, y = value, color = chain)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  labs(title = "MCMC Trace Plots (1-Pool Model)",
       subtitle = "Each color is a separate chain. Look for good mixing.",
       x = "Iteration", y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(trace_1pool)

# =============================================================================
# STEP 1.7: Remove burn-in and thin the chains
# =============================================================================
# BURN-IN: The early iterations before chains have found the high-probability
#   region. We discard these.
#
# THINNING: Keep every Nth sample to reduce autocorrelation (successive
#   samples being too similar). This gives more independent draws.
#
# *** EXERCISE: Try different burn-in and thinning values! ***
#   - What happens with no burn-in (burn_in = 1)?
#   - What if you thin too aggressively (thin = 100)?
#   - What if you don't thin at all (thin = 1)?

burn_in_1pool <- 2000   # discard first 2000 iterations
thin_1pool    <- 5      # keep every 5th sample

samples_1pool <- getSample(mcmc_1pool, start = burn_in_1pool, thin = thin_1pool)
colnames(samples_1pool) <- c("k", "input")

cat("\nRetained", nrow(samples_1pool), "posterior samples after burn-in & thinning\n")

# Visualize burn-in cutoff on trace plot
trace_burnin_1 <- ggplot(chain_long_1, aes(x = iteration, y = value, color = chain)) +
  geom_line(alpha = 0.5) +
  geom_vline(xintercept = burn_in_1pool / 3, color = "black",
             linewidth = 1, linetype = "dotted") +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  labs(title = "Burn-in Removal",
       subtitle = "Everything left of the dotted line is discarded",
       x = "Iteration", y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(trace_burnin_1)

# =============================================================================
# STEP 1.8: Compare prior and posterior
# =============================================================================
# The posterior is the prior updated by the data. If the data are informative,
# the posterior will be much narrower than the prior.

posterior_df_1 <- as.data.frame(samples_1pool)

# Posterior summaries
post_summary_1 <- posterior_df_1 %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    median   = median(value),
    lower_90 = quantile(value, 0.05),
    upper_90 = quantile(value, 0.95),
    .groups  = "drop"
  ) %>%
  mutate(
    prior_lo = prior_lower_1pool[parameter],
    prior_hi = prior_upper_1pool[parameter],
    uncert_reduction = round(100 * (1 - (upper_90 - lower_90) / (prior_hi - prior_lo)), 1)
  )

map_1pool <- MAP(mcmc_1pool, start = burn_in_1pool)$parametersMAP
names(map_1pool) <- c("k", "input")

cat("\n=== 1-Pool Posterior Summary ===\n")
cat(sprintf("%-8s  %8s  %10s  %10s  %10s\n",
            "Param", "MAP", "Median", "90% lo", "90% hi"))
for (i in seq_len(nrow(post_summary_1))) {
  r <- post_summary_1[i, ]
  cat(sprintf("%-8s  %8.4f  %10.4f  %10.4f  %10.4f  (%.0f%% reduction)\n",
              r$parameter, map_1pool[r$parameter], r$median,
              r$lower_90, r$upper_90, r$uncert_reduction))
}

# Plot prior vs posterior
posterior_long_1 <- posterior_df_1 %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c("k", "input"),
                            labels = c("k (yr^-1)", "input (gC/m2/yr)")))

prior_vs_post_1 <- ggplot() +
  geom_histogram(data = prior_long_1pool,
                 aes(x = value, y = after_stat(density)),
                 bins = 50, fill = "skyblue", alpha = 0.4, color = "white") +
  geom_density(data = posterior_long_1, aes(x = value),
               fill = "coral", alpha = 0.5, linewidth = 1) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +
  labs(title = "Prior (blue) vs. Posterior (coral) --1-Pool Model",
       subtitle = "How much did the data narrow down our parameter estimates?",
       x = "Parameter value", y = "Density") +
  theme_minimal(base_size = 14)

print(prior_vs_post_1)

# =============================================================================
# STEP 1.9: Plot model predictions using the posterior
# =============================================================================
# Run the model with the MAP estimate and an ensemble of posterior draws
# to show both the best fit and the uncertainty.

years_plot_1 <- seq(2009, 2022, by = 0.5)

# MAP prediction
map_model_1 <- OnepModel14(
  t = years_plot_1, k = map_1pool["k"], C0 = C0_1pool,
  F0_Delta14C = F0_1pool, In = map_1pool["input"],
  inputFc = atmIn, pass = TRUE
)

map_pred_1 <- data.frame(
  year    = years_plot_1,
  C14     = as.numeric(getF14(map_model_1)),
  C_stock = as.numeric(getC(map_model_1))
)

# Posterior ensemble for uncertainty envelope
n_draws_1 <- min(200, nrow(posterior_df_1))
draw_idx_1 <- sample(nrow(posterior_df_1), n_draws_1)

ensemble_1 <- bind_rows(lapply(draw_idx_1, function(i) {
  m <- tryCatch(
    OnepModel14(
      t = years_plot_1, k = posterior_df_1$k[i], C0 = C0_1pool,
      F0_Delta14C = F0_1pool, In = posterior_df_1$input[i],
      inputFc = atmIn, pass = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  data.frame(year = years_plot_1,
             C14 = as.numeric(getF14(m)),
             C_stock = as.numeric(getC(m)))
}))

envelope_1 <- ensemble_1 %>%
  group_by(year) %>%
  summarise(
    C14_lo = quantile(C14, 0.05), C14_hi = quantile(C14, 0.95),
    Cstock_lo = quantile(C_stock, 0.05), Cstock_hi = quantile(C_stock, 0.95),
    .groups = "drop"
  )

# Plot
p1_c14 <- ggplot() +
  geom_ribbon(data = envelope_1, aes(x = year, ymin = C14_lo, ymax = C14_hi),
              fill = "steelblue", alpha = 0.3) +
  geom_line(data = map_pred_1, aes(x = year, y = C14),
            color = "steelblue", linewidth = 1.2) +
  geom_point(data = dat_1pool, aes(x = year, y = C14), size = 4) +
  geom_errorbar(data = dat_1pool,
                aes(x = year, ymin = C14 - C_14_sd, ymax = C14 + C_14_sd),
                width = 0.3, linewidth = 0.8) +
  labs(title = expression(paste("Bulk Soil ", Delta^14, "C")),
       subtitle = "Line = MAP estimate. Shaded = 90% posterior envelope.",
       x = "Year", y = expression(paste(Delta^14, "C (per mil)"))) +
  theme_minimal(base_size = 14)

p1_cstock <- ggplot() +
  geom_ribbon(data = envelope_1, aes(x = year, ymin = Cstock_lo/1000, ymax = Cstock_hi/1000),
              fill = "forestgreen", alpha = 0.3) +
  geom_line(data = map_pred_1, aes(x = year, y = C_stock/1000),
            color = "forestgreen", linewidth = 1.2) +
  geom_point(data = dat_1pool, aes(x = year, y = C_stock/1000), size = 4) +
  geom_errorbar(data = dat_1pool,
                aes(x = year, ymin = (C_stock - C_stock_sd)/1000,
                    ymax = (C_stock + C_stock_sd)/1000),
                width = 0.3, linewidth = 0.8) +
  labs(title = "Bulk Soil Carbon Stock",
       subtitle = "Line = MAP estimate. Shaded = 90% posterior envelope.",
       x = "Year", y = expression(paste("C stock (kg m"^-2, ")"))) +
  theme_minimal(base_size = 14)

print(p1_c14 | p1_cstock)

cat("\n*** 1-Pool model complete! ***\n")
cat("Before moving on, try the exercises:\n")
cat("  - Change the prior bounds (STEP 1.3) and re-run from there\n")
cat("  - Change burn-in and thinning (STEP 1.7)\n")
cat("  - Change the number of MCMC iterations (STEP 1.5)\n\n")


# #############################################################################
#
#  PART 2: THREE-POOL SERIES MODEL  (6 parameters)
#
# #############################################################################
#
# Now we model the three soil layers explicitly as pools connected in series:
#   Surface organic -> Deep organic -> Mineral
#
# Carbon enters the surface pool, decomposes in each pool, and transfers
# downward through the profile. This is more realistic but much harder
# to fit --6 parameters can trade off against each other.
#
# Parameters to estimate:
#   k1, k2, k3  = decomposition rates for surface, deep, mineral (yr^-1)
#   input        = total carbon input rate to surface pool (gC/m2/yr)
#   a21          = transfer rate from surface (pool 1) to deep (pool 2)
#   a32          = transfer rate from deep (pool 2) to mineral (pool 3)

# =============================================================================
# STEP 2.1: Prepare 3-pool observations
# =============================================================================
# Each layer is now its own pool. Use Control treatment.

dat_3pool <- datComb %>%
  filter(year == 2009 | treatment == "Control") %>%
  mutate(C_stock    = C_stock * 1000,
         C_stock_sd = C_stock_sd * 1000)

# Initial conditions from 2009 (ordered: so, do, min)
init_3pool  <- dat_3pool %>% filter(year == 2009)
final_3pool <- dat_3pool %>% filter(year == 2022)

C0_3pool <- init_3pool$C_stock
F0_3pool <- init_3pool$C14
years_3pool <- seq(2009, 2022)

cat("\n=== 3-Pool Initial Conditions (2009) ===\n")
cat("Layers:", init_3pool$layer, "\n")
cat("C stocks (gC/m2):", round(C0_3pool), "\n")
cat("Delta14C (per mil):", round(F0_3pool, 1), "\n")

# =============================================================================
# STEP 2.2: Define the likelihood function
# =============================================================================
# Same idea as the 1-pool model, but now we compare predictions for each
# of the 3 pools separately --that's 6 data points (3 C stocks + 3 14C values)
# constraining 6 parameters.

likelihood_3pool <- function(pars) {
  k1 <- pars[1]; k2 <- pars[2]; k3 <- pars[3]
  input <- pars[4]; a21 <- pars[5]; a32 <- pars[6]

  model <- tryCatch(
    ThreepSeriesModel14(
      t = years_3pool,
      ks = c(k1 = k1, k2 = k2, k3 = k3),
      C0 = C0_3pool,
      F0_Delta14C = F0_3pool,
      In = input,
      a21 = a21,
      a32 = a32,
      inputFc = atmIn,
      pass = TRUE
    ),
    error = function(e) return(NULL)
  )

  if (is.null(model)) return(-1e10)

  idx <- length(years_3pool)
  pred_C14    <- getF14(model)[idx, ]
  pred_Cstock <- getC(model)[idx, ]

  # Compare each pool to its observed value
  ll_C14    <- sum(dnorm(pred_C14 - final_3pool$C14, 0,
                         final_3pool$C_14_sd * 0.4, log = TRUE))
  ll_Cstock <- sum(dnorm(pred_Cstock - final_3pool$C_stock, 0,
                         final_3pool$C_stock_sd * 0.4, log = TRUE))

  return(sum(ll_C14, ll_Cstock, na.rm = TRUE))
}

# Sanity check
cat("\nTest likelihood at example parameters:",
    likelihood_3pool(c(0.1, 0.01, 0.01, 200, 0.1, 0.01)), "\n")

# =============================================================================
# STEP 2.3: Define priors
# =============================================================================
# With 6 parameters, prior choice matters a lot! Too wide and the MCMC
# wastes time in impossible regions. Too narrow and you might miss the answer.
#
# These ranges are informed by the literature on permafrost soils:
#   - Surface organic decomposes fastest (k1 largest)
#   - Mineral soil decomposes slowest (k3 smallest)
#   - Transfer rates must be smaller than decomposition rates
#
# *** EXERCISE: Try different priors! ***
#   - What if you swap k1 and k3 ranges (allow slow surface, fast mineral)?
#   - What happens with much wider priors?
#   - What if you constrain input very tightly?

prior_lower_3pool <- c(k1 = 0.002,  k2 = 0.0006, k3 = 0.0002,
                        input = 10,  a21 = 0.0025, a32 = 0.0001)
prior_upper_3pool <- c(k1 = 0.2,    k2 = 0.06,   k3 = 0.06,
                        input = 205, a21 = 0.5,    a32 = 0.1)

prior_3pool <- createUniformPrior(
  lower = prior_lower_3pool,
  upper = prior_upper_3pool
)

# =============================================================================
# STEP 2.4: Visualize the 3-pool priors
# =============================================================================

par_names_3pool <- names(prior_lower_3pool)

prior_samples_3pool <- as.data.frame(mapply(
  runif, n = 5000, min = prior_lower_3pool, max = prior_upper_3pool
))

prior_long_3pool <- prior_samples_3pool %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = par_names_3pool))

prior_3pool_plot <- ggplot(prior_long_3pool, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "skyblue", alpha = 0.7, color = "white") +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(title = "Prior Distributions (3-Pool Model)",
       x = "Parameter value", y = "Density") +
  theme_minimal(base_size = 12)

print(prior_3pool_plot)

# =============================================================================
# STEP 2.5: Run the MCMC
# =============================================================================
# 6 parameters need more iterations. With only 13 years of data and 6 data
# points constraining 6 parameters, this is a hard problem --expect wider
# posteriors and more parameter correlations than the 1-pool model.

setup_3pool <- createBayesianSetup(
  likelihood = likelihood_3pool,
  prior      = prior_3pool,
  names      = par_names_3pool
)

# 20,000 iterations for a teaching example. For publication you'd want more
# (the original analysis used 500,000 for some treatments).
settings_3pool <- list(iterations = 20000, message = FALSE)

cat("\nRunning 3-pool MCMC (this may take 1-5 minutes)...\n")
mcmc_3pool <- runMCMC(
  bayesianSetup = setup_3pool,
  sampler       = "DREAMzs",
  settings      = settings_3pool
)
cat("Done!\n")

# =============================================================================
# STEP 2.6: Trace plots for all 6 parameters
# =============================================================================

raw_chains_3 <- getSample(mcmc_3pool, start = 1, coda = TRUE)

chain_df_3 <- do.call(rbind, lapply(seq_along(raw_chains_3), function(i) {
  ch <- as.data.frame(as.matrix(raw_chains_3[[i]]))
  colnames(ch) <- par_names_3pool
  ch$iteration <- seq_len(nrow(ch))
  ch$chain <- factor(i)
  ch
}))

chain_long_3 <- chain_df_3 %>%
  pivot_longer(cols = all_of(par_names_3pool),
               names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = par_names_3pool))

trace_3pool <- ggplot(chain_long_3, aes(x = iteration, y = value, color = chain)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(title = "MCMC Trace Plots (3-Pool Model)",
       subtitle = "Look for mixing and convergence across all chains",
       x = "Iteration", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

print(trace_3pool)

# =============================================================================
# STEP 2.7: Remove burn-in and thin
# =============================================================================
# With more parameters, you typically need a longer burn-in.
#
# *** EXERCISE: Look at your trace plots and adjust these! ***
#   - If chains haven't converged by 5000, increase burn_in
#   - If you see high autocorrelation, increase thin

burn_in_3pool <- 5000
thin_3pool    <- 10

samples_3pool <- getSample(mcmc_3pool, start = burn_in_3pool, thin = thin_3pool)
colnames(samples_3pool) <- par_names_3pool
samples_df_3 <- as.data.frame(samples_3pool)

cat("\nRetained", nrow(samples_df_3), "posterior samples after burn-in & thinning\n")

# =============================================================================
# STEP 2.8: Prior vs. Posterior for all parameters
# =============================================================================

posterior_long_3 <- samples_df_3 %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = par_names_3pool))

# Posterior summaries
post_summary_3 <- posterior_long_3 %>%
  group_by(parameter) %>%
  summarise(
    median   = median(value),
    lower_90 = quantile(value, 0.05),
    upper_90 = quantile(value, 0.95),
    .groups  = "drop"
  ) %>%
  mutate(
    prior_range = prior_upper_3pool[parameter] - prior_lower_3pool[parameter],
    post_range  = upper_90 - lower_90,
    uncert_reduction = round(100 * (1 - post_range / prior_range), 1)
  )

map_3pool <- MAP(mcmc_3pool, start = burn_in_3pool)$parametersMAP
names(map_3pool) <- par_names_3pool

cat("\n=== 3-Pool Posterior Summary ===\n")
cat(sprintf("%-8s  %8s  %10s  %10s  %10s  %s\n",
            "Param", "MAP", "Median", "90% lo", "90% hi", "Reduction"))
for (i in seq_len(nrow(post_summary_3))) {
  r <- post_summary_3[i, ]
  cat(sprintf("%-8s  %8.4f  %10.4f  %10.4f  %10.4f  %5.1f%%\n",
              r$parameter, map_3pool[r$parameter], r$median,
              r$lower_90, r$upper_90, r$uncert_reduction))
}

# Combined prior + posterior density plot
prior_vs_post_3 <- ggplot() +
  geom_histogram(data = prior_long_3pool,
                 aes(x = value, y = after_stat(density)),
                 bins = 40, fill = "skyblue", alpha = 0.4, color = "white") +
  geom_density(data = posterior_long_3, aes(x = value),
               fill = "coral", alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = data.frame(
    parameter = factor(par_names_3pool, levels = par_names_3pool),
    xval = map_3pool),
    aes(xintercept = xval), linetype = "dashed", color = "darkblue") +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(title = "Prior (blue) vs. Posterior (coral) --3-Pool Model",
       subtitle = "Dashed blue = MAP estimate. Narrow posteriors = well-constrained parameters.",
       x = "Parameter value", y = "Density") +
  theme_minimal(base_size = 12)

print(prior_vs_post_3)

# =============================================================================
# STEP 2.9: Posterior correlations
# =============================================================================
# In multi-parameter models, parameters often trade off against each other.
# Strong correlations mean the data can't independently constrain those
# parameters --only certain combinations are well determined.

cat("\n--- Posterior Correlation Matrix ---\n")
print(round(cor(samples_df_3), 2))

# Pairs plot (subsample for speed)
n_pairs <- min(500, nrow(samples_df_3))
pairs(samples_df_3[sample(nrow(samples_df_3), n_pairs), ],
      pch = 16, cex = 0.4, col = adjustcolor("steelblue", 0.3),
      main = "Posterior Pairwise Correlations (3-Pool Model)")

# =============================================================================
# STEP 2.10: Plot model predictions with posterior uncertainty
# =============================================================================
# Forward-run the model at the MAP and with posterior draws to visualize
# how well the calibrated model fits the observations.

pool_labels <- c(so = "Surface organic", do = "Deep organic", min = "Mineral")
years_plot_3 <- seq(2009, 2022, by = 0.5)

# MAP forward run
map_model_3 <- ThreepSeriesModel14(
  t = years_plot_3,
  ks = c(k1 = map_3pool["k1"], k2 = map_3pool["k2"], k3 = map_3pool["k3"]),
  C0 = C0_3pool, F0_Delta14C = F0_3pool,
  In = map_3pool["input"], a21 = map_3pool["a21"], a32 = map_3pool["a32"],
  inputFc = atmIn, pass = TRUE
)

map_pred_3 <- data.frame(
  year = rep(years_plot_3, 3),
  pool = rep(pool_labels, each = length(years_plot_3)),
  C14     = as.vector(getF14(map_model_3)),
  C_stock = as.vector(getC(map_model_3))
)

# Posterior ensemble
n_draws_3 <- min(200, nrow(samples_df_3))
draw_idx_3 <- sample(nrow(samples_df_3), n_draws_3)

ensemble_3 <- bind_rows(lapply(draw_idx_3, function(i) {
  p <- as.numeric(samples_df_3[i, ])
  m <- tryCatch(
    ThreepSeriesModel14(
      t = years_plot_3,
      ks = c(k1 = p[1], k2 = p[2], k3 = p[3]),
      C0 = C0_3pool, F0_Delta14C = F0_3pool,
      In = p[4], a21 = p[5], a32 = p[6],
      inputFc = atmIn, pass = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  data.frame(
    year = rep(years_plot_3, 3),
    pool = rep(pool_labels, each = length(years_plot_3)),
    C14     = as.vector(getF14(m)),
    C_stock = as.vector(getC(m))
  )
}))

envelope_3 <- ensemble_3 %>%
  group_by(year, pool) %>%
  summarise(
    C14_lo = quantile(C14, 0.05), C14_hi = quantile(C14, 0.95),
    Cstock_lo = quantile(C_stock, 0.05), Cstock_hi = quantile(C_stock, 0.95),
    .groups = "drop"
  )

# Observation data for plotting
obs_plot_3 <- dat_3pool %>%
  mutate(pool = factor(pool_labels[layer], levels = pool_labels))

map_pred_3$pool  <- factor(map_pred_3$pool, levels = pool_labels)
envelope_3$pool  <- factor(envelope_3$pool, levels = pool_labels)

# 14C by pool
p3_c14 <- ggplot() +
  geom_ribbon(data = envelope_3,
              aes(x = year, ymin = C14_lo, ymax = C14_hi),
              fill = "steelblue", alpha = 0.3) +
  geom_line(data = map_pred_3, aes(x = year, y = C14),
            color = "steelblue", linewidth = 1) +
  geom_point(data = obs_plot_3, aes(x = year, y = C14), size = 3) +
  geom_errorbar(data = obs_plot_3,
                aes(x = year, ymin = C14 - C_14_sd, ymax = C14 + C_14_sd),
                width = 0.3) +
  facet_wrap(~ pool, ncol = 1, scales = "free_y") +
  labs(title = expression(paste(Delta^14, "C by Soil Layer")),
       subtitle = "MAP prediction + 90% posterior envelope vs. observations",
       x = "Year", y = expression(paste(Delta^14, "C (per mil)"))) +
  theme_minimal(base_size = 12)

# C stock by pool
p3_cstock <- ggplot() +
  geom_ribbon(data = envelope_3,
              aes(x = year, ymin = Cstock_lo/1000, ymax = Cstock_hi/1000),
              fill = "forestgreen", alpha = 0.3) +
  geom_line(data = map_pred_3, aes(x = year, y = C_stock/1000),
            color = "forestgreen", linewidth = 1) +
  geom_point(data = obs_plot_3, aes(x = year, y = C_stock/1000), size = 3) +
  geom_errorbar(data = obs_plot_3,
                aes(x = year, ymin = (C_stock - C_stock_sd)/1000,
                    ymax = (C_stock + C_stock_sd)/1000),
                width = 0.3) +
  facet_wrap(~ pool, ncol = 1, scales = "free_y") +
  labs(title = "Carbon Stock by Soil Layer",
       subtitle = "MAP prediction + 90% posterior envelope vs. observations",
       x = "Year", y = expression(paste("C stock (kg m"^-2, ")"))) +
  theme_minimal(base_size = 12)

print(p3_c14 | p3_cstock)

# =============================================================================
# STEP 2.11: Summary
# =============================================================================

cat("\n=== Uncertainty Reduction Summary (3-Pool) ===\n")
cat(sprintf("%-8s  %8s  %10s  %10s  %s\n",
            "Param", "MAP", "90% lo", "90% hi", "Reduction"))
for (i in seq_len(nrow(post_summary_3))) {
  r <- post_summary_3[i, ]
  cat(sprintf("%-8s  %8.4f  %10.4f  %10.4f  %5.1f%%\n",
              r$parameter, map_3pool[r$parameter],
              r$lower_90, r$upper_90, r$uncert_reduction))
}

cat("\n=== Tutorial complete! ===\n")
cat("Things to try:\n")
cat("  1. Change the prior bounds and see how posteriors shift\n")
cat("  2. Increase/decrease MCMC iterations and watch convergence\n")
cat("  3. Change burn-in and thinning --how does it affect results?\n")
cat("  4. Scale the observation uncertainty (the 0.4 factor in the likelihood)\n")
cat("  5. Try the Warming treatment instead of Control\n")
