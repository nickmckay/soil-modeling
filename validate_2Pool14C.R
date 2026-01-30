# =============================================================================
# Validation: Two-Pool Model - SoilR vs N-Pool Implementation
# =============================================================================
#
# This script validates the N-pool model implementation against SoilR's
# two-pool models by running identical experiments and comparing results.
#
# We compare three model structures using both implementations:
#   - Parallel: Independent pools with direct atmospheric input
#   - Series: Sequential flow from fast to slow pool
#   - Feedback: Series with backward transfer
#
# Parameters match exactly those in 2Pool14C.R
# =============================================================================

library(SoilR)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# PART 1: Load N-Pool Model Functions
# =============================================================================

source("NPool14C.R")

# =============================================================================
# PART 2: Set Up Common Parameters
# =============================================================================
# These parameters match exactly those in 2Pool14C.R

# Atmospheric 14C data
ad <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone1, time.scale = "AD")
Fatm <- ad[, 1:2]
atm_data <- prepare_atm_data(Fatm)

# Time vector
years <- seq(1901, 2009, by = 0.5)

# Annual litter input
LitterInput <- 100

# Decomposition rate constants
# Pool 1 (fast): k1 = 1/2 yr^-1 → turnover time = 2 years
# Pool 2 (slow): k2 = 1/50 yr^-1 → turnover time = 50 years
k <- c(k1 = 1/2, k2 = 1/50)

# Initial carbon stocks
C0 <- c(100, 1000)

# Initial 14C (pre-bomb = 0 per mil)
F0_Delta14C <- c(0, 0)

# Transfer coefficients
gam <- 0.7           # 70% of input to fast pool (parallel model)
a21 <- 0.5 * k[1]    # 50% transfer from pool 1 to pool 2
a12 <- 0.3 * k[2]    # 30% feedback from pool 2 to pool 1

# =============================================================================
# PART 3: Run SoilR Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Running SoilR Two-Pool Models...\n")
cat("=======================================================================\n")

# Parallel model
soilr_parallel <- TwopParallelModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  gam = gam,
  inputFc = Fatm
)

# Series model
soilr_series <- TwopSeriesModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  a21 = a21,
  inputFc = Fatm
)

# Feedback model
soilr_feedback <- TwopFeedbackModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  a21 = a21,
  a12 = a12,
  inputFc = Fatm
)

# Extract results
soilr_par_pools <- getF14(soilr_parallel)
soilr_par_bulk <- getF14C(soilr_parallel)
soilr_par_resp <- getF14R(soilr_parallel)

soilr_ser_pools <- getF14(soilr_series)
soilr_ser_bulk <- getF14C(soilr_series)
soilr_ser_resp <- getF14R(soilr_series)

soilr_fb_pools <- getF14(soilr_feedback)
soilr_fb_bulk <- getF14C(soilr_feedback)
soilr_fb_resp <- getF14R(soilr_feedback)

cat("  SoilR models complete.\n")

# =============================================================================
# PART 4: Run N-Pool Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Running N-Pool Two-Pool Models...\n")
cat("=======================================================================\n")

# Parallel model: no transfers between pools
transfer_parallel <- create_parallel_transfers(2)

npool_parallel <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(gam, 1 - gam),
  transfer_matrix = transfer_parallel,
  F_atm_data = atm_data
)

# Series model: transfer from pool 1 to pool 2, all input to pool 1
transfer_series <- matrix(0, 2, 2)
transfer_series[2, 1] <- a21

npool_series <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(1, 0),  # All input to pool 1
  transfer_matrix = transfer_series,
  F_atm_data = atm_data
)

# Feedback model: bidirectional transfers
transfer_feedback <- matrix(0, 2, 2)
transfer_feedback[2, 1] <- a21  # Forward: pool 1 -> pool 2
transfer_feedback[1, 2] <- a12  # Backward: pool 2 -> pool 1

npool_feedback <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(1, 0),
  transfer_matrix = transfer_feedback,
  F_atm_data = atm_data
)

# Extract results
npool_par_pools <- getF14_npool(npool_parallel)
npool_par_bulk <- getF14C_npool(npool_parallel)
npool_par_resp <- getF14R_npool(npool_parallel)

npool_ser_pools <- getF14_npool(npool_series)
npool_ser_bulk <- getF14C_npool(npool_series)
npool_ser_resp <- getF14R_npool(npool_series)

npool_fb_pools <- getF14_npool(npool_feedback)
npool_fb_bulk <- getF14C_npool(npool_feedback)
npool_fb_resp <- getF14R_npool(npool_feedback)

cat("  N-Pool models complete.\n")

# =============================================================================
# PART 5: Calculate Differences
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Validation Results: Two-Pool Models\n")
cat("=======================================================================\n")

# Parallel model
cat("\nPARALLEL MODEL:\n")
cat(sprintf("  Pool 1 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_par_pools[,1] - npool_par_pools[,1])),
            sqrt(mean((soilr_par_pools[,1] - npool_par_pools[,1])^2))))
cat(sprintf("  Pool 2 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_par_pools[,2] - npool_par_pools[,2])),
            sqrt(mean((soilr_par_pools[,2] - npool_par_pools[,2])^2))))
cat(sprintf("  Bulk SOM - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_par_bulk - npool_par_bulk)),
            sqrt(mean((soilr_par_bulk - npool_par_bulk)^2))))
cat(sprintf("  Respired - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_par_resp - npool_par_resp)),
            sqrt(mean((soilr_par_resp - npool_par_resp)^2))))

# Series model
cat("\nSERIES MODEL:\n")
cat(sprintf("  Pool 1 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_ser_pools[,1] - npool_ser_pools[,1])),
            sqrt(mean((soilr_ser_pools[,1] - npool_ser_pools[,1])^2))))
cat(sprintf("  Pool 2 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_ser_pools[,2] - npool_ser_pools[,2])),
            sqrt(mean((soilr_ser_pools[,2] - npool_ser_pools[,2])^2))))
cat(sprintf("  Bulk SOM - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_ser_bulk - npool_ser_bulk)),
            sqrt(mean((soilr_ser_bulk - npool_ser_bulk)^2))))
cat(sprintf("  Respired - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_ser_resp - npool_ser_resp)),
            sqrt(mean((soilr_ser_resp - npool_ser_resp)^2))))

# Feedback model
cat("\nFEEDBACK MODEL:\n")
cat(sprintf("  Pool 1 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_fb_pools[,1] - npool_fb_pools[,1])),
            sqrt(mean((soilr_fb_pools[,1] - npool_fb_pools[,1])^2))))
cat(sprintf("  Pool 2 - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_fb_pools[,2] - npool_fb_pools[,2])),
            sqrt(mean((soilr_fb_pools[,2] - npool_fb_pools[,2])^2))))
cat(sprintf("  Bulk SOM - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_fb_bulk - npool_fb_bulk)),
            sqrt(mean((soilr_fb_bulk - npool_fb_bulk)^2))))
cat(sprintf("  Respired - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
            max(abs(soilr_fb_resp - npool_fb_resp)),
            sqrt(mean((soilr_fb_resp - npool_fb_resp)^2))))

# =============================================================================
# PART 6: Prepare Data for Visualization
# =============================================================================

# Atmospheric data
atm_df <- data.frame(Year = Fatm[, 1], Delta14C = Fatm[, 2])

# Helper function for pool data
create_pool_comparison_df <- function(years, soilr_pools, npool_pools, model_name) {
  bind_rows(
    data.frame(Year = years, Delta14C = soilr_pools[, 1],
               Pool = "Pool 1 (fast)", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = soilr_pools[, 2],
               Pool = "Pool 2 (slow)", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = npool_pools[, 1],
               Pool = "Pool 1 (fast)", Source = "N-Pool", Model = model_name),
    data.frame(Year = years, Delta14C = npool_pools[, 2],
               Pool = "Pool 2 (slow)", Source = "N-Pool", Model = model_name)
  )
}

# Helper function for bulk/respired data
create_bulk_comparison_df <- function(years, soilr_bulk, soilr_resp,
                                       npool_bulk, npool_resp, model_name) {
  bind_rows(
    data.frame(Year = years, Delta14C = soilr_bulk,
               Type = "Bulk SOM", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = soilr_resp,
               Type = "Respired CO2", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = npool_bulk,
               Type = "Bulk SOM", Source = "N-Pool", Model = model_name),
    data.frame(Year = years, Delta14C = npool_resp,
               Type = "Respired CO2", Source = "N-Pool", Model = model_name)
  )
}

# Combine pool data
pools_df <- bind_rows(
  create_pool_comparison_df(years, soilr_par_pools, npool_par_pools, "Parallel"),
  create_pool_comparison_df(years, soilr_ser_pools, npool_ser_pools, "Series"),
  create_pool_comparison_df(years, soilr_fb_pools, npool_fb_pools, "Feedback")
)
pools_df$Model <- factor(pools_df$Model, levels = c("Parallel", "Series", "Feedback"))

# Combine bulk/respired data
bulk_df <- bind_rows(
  create_bulk_comparison_df(years, soilr_par_bulk, soilr_par_resp,
                            npool_par_bulk, npool_par_resp, "Parallel"),
  create_bulk_comparison_df(years, soilr_ser_bulk, soilr_ser_resp,
                            npool_ser_bulk, npool_ser_resp, "Series"),
  create_bulk_comparison_df(years, soilr_fb_bulk, soilr_fb_resp,
                            npool_fb_bulk, npool_fb_resp, "Feedback")
)
bulk_df$Model <- factor(bulk_df$Model, levels = c("Parallel", "Series", "Feedback"))

# Difference data for bulk SOM
diff_bulk_df <- bind_rows(
  data.frame(Year = years, Difference = soilr_par_bulk - npool_par_bulk, Model = "Parallel"),
  data.frame(Year = years, Difference = soilr_ser_bulk - npool_ser_bulk, Model = "Series"),
  data.frame(Year = years, Difference = soilr_fb_bulk - npool_fb_bulk, Model = "Feedback")
)
diff_bulk_df$Model <- factor(diff_bulk_df$Model, levels = c("Parallel", "Series", "Feedback"))

# =============================================================================
# PART 7: Visualization
# =============================================================================

# Plot 1: Individual pools comparison
p_pools <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(pools_df, Source == "SoilR"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "black", linewidth = 0.9) +
  geom_line(data = filter(pools_df, Source == "N-Pool"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "red", linewidth = 0.7) +
  facet_wrap(~Model, ncol = 1) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Individual Pool Comparison",
    subtitle = "Black = SoilR | Red = N-Pool",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Plot 2: Bulk SOM and Respired comparison
p_bulk <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(bulk_df, Source == "SoilR"),
            aes(x = Year, y = Delta14C, color = Type),
            linewidth = 0.9) +
  geom_line(data = filter(bulk_df, Source == "N-Pool"),
            aes(x = Year, y = Delta14C, color = Type),
            linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~Model, ncol = 1) +
  scale_color_manual(values = c("Bulk SOM" = "blue", "Respired CO2" = "red")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Bulk SOM & Respired CO2",
    subtitle = "Solid = SoilR | Dashed = N-Pool",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    color = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Plot 3: Difference plot for bulk SOM
p_diff <- ggplot(diff_bulk_df, aes(x = Year, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "darkblue", linewidth = 0.8) +
  facet_wrap(~Model, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Bulk SOM Difference",
    subtitle = "SoilR - N-Pool",
    x = "Year",
    y = expression(Delta*Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"))

# Display plots
print(p_pools)
print(p_bulk)
print(p_diff)

# Combined figure
(p_pools | p_bulk | p_diff) +
  plot_annotation(
    title = "Two-Pool Model Validation: SoilR vs N-Pool",
    subtitle = "Comparing Parallel, Series, and Feedback structures"
  )
