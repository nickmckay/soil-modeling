# =============================================================================
# Validation: Three-Pool Model - SoilR vs N-Pool Implementation
# =============================================================================
#
# This script validates the N-pool model implementation against SoilR's
# three-pool models by running identical experiments and comparing results.
#
# We compare three model structures using both implementations:
#   - Parallel: Independent pools with direct atmospheric input
#   - Series: Sequential flow fast → intermediate → slow
#   - Feedback: Series with backward transfers
#
# Parameters match exactly those in 3Pool14C.R
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
# These parameters match exactly those in 3Pool14C.R

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
# Pool 2 (intermediate): k2 = 1/10 yr^-1 → turnover time = 10 years
# Pool 3 (slow): k3 = 1/50 yr^-1 → turnover time = 50 years
k <- c(k1 = 1/2, k2 = 1/10, k3 = 1/50)

# Initial carbon stocks
C0 <- c(100, 500, 1000)

# Initial 14C (pre-bomb = 0 per mil)
F0_Delta14C <- c(0, 0, 0)

# Input allocation for parallel model
gam1 <- 0.6  # 60% to fast pool
gam2 <- 0.2  # 20% to intermediate pool
# (20% to slow pool)

# Transfer coefficients for series/feedback models
a21 <- 0.9 * k[1]  # 90% transfer from pool 1 to pool 2
a32 <- 0.4 * k[2]  # 40% transfer from pool 2 to pool 3
a12 <- 0.4 * k[2]  # 40% feedback from pool 2 to pool 1
a23 <- 0.7 * k[3]  # 70% feedback from pool 3 to pool 2

# =============================================================================
# PART 3: Run SoilR Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Running SoilR Three-Pool Models...\n")
cat("=======================================================================\n")

# Parallel model
soilr_parallel <- ThreepParallelModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  gam1 = gam1,
  gam2 = gam2,
  inputFc = Fatm
)

# Series model
soilr_series <- ThreepSeriesModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  a21 = a21,
  a32 = a32,
  inputFc = Fatm
)

# Feedback model
soilr_feedback <- ThreepFeedbackModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  a21 = a21,
  a12 = a12,
  a32 = a32,
  a23 = a23,
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
cat("Running N-Pool Three-Pool Models...\n")
cat("=======================================================================\n")

# Parallel model: no transfers between pools
transfer_parallel <- create_parallel_transfers(3)

npool_parallel <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(gam1, gam2, 1 - gam1 - gam2),
  transfer_matrix = transfer_parallel,
  F_atm_data = atm_data
)

# Series model: sequential transfers, all input to pool 1
transfer_series <- matrix(0, 3, 3)
transfer_series[2, 1] <- a21  # Pool 1 -> Pool 2
transfer_series[3, 2] <- a32  # Pool 2 -> Pool 3

npool_series <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(1, 0, 0),  # All input to pool 1
  transfer_matrix = transfer_series,
  F_atm_data = atm_data
)

# Feedback model: bidirectional transfers
transfer_feedback <- matrix(0, 3, 3)
transfer_feedback[2, 1] <- a21  # Forward: Pool 1 -> Pool 2
transfer_feedback[1, 2] <- a12  # Backward: Pool 2 -> Pool 1
transfer_feedback[3, 2] <- a32  # Forward: Pool 2 -> Pool 3
transfer_feedback[2, 3] <- a23  # Backward: Pool 3 -> Pool 2

npool_feedback <- solve_npool_model(
  times = years,
  k = k,
  C0 = C0,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = c(1, 0, 0),
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
cat("Validation Results: Three-Pool Models\n")
cat("=======================================================================\n")

# Function to print validation stats
print_validation_stats <- function(name, soilr_val, npool_val) {
  diff <- soilr_val - npool_val
  cat(sprintf("  %s - Max diff: %.6f ‰, RMSE: %.6f ‰\n",
              name, max(abs(diff)), sqrt(mean(diff^2))))
}

# Parallel model
cat("\nPARALLEL MODEL:\n")
print_validation_stats("Pool 1 (fast)", soilr_par_pools[,1], npool_par_pools[,1])
print_validation_stats("Pool 2 (intermediate)", soilr_par_pools[,2], npool_par_pools[,2])
print_validation_stats("Pool 3 (slow)", soilr_par_pools[,3], npool_par_pools[,3])
print_validation_stats("Bulk SOM", soilr_par_bulk, npool_par_bulk)
print_validation_stats("Respired CO2", soilr_par_resp, npool_par_resp)

# Series model
cat("\nSERIES MODEL:\n")
print_validation_stats("Pool 1 (fast)", soilr_ser_pools[,1], npool_ser_pools[,1])
print_validation_stats("Pool 2 (intermediate)", soilr_ser_pools[,2], npool_ser_pools[,2])
print_validation_stats("Pool 3 (slow)", soilr_ser_pools[,3], npool_ser_pools[,3])
print_validation_stats("Bulk SOM", soilr_ser_bulk, npool_ser_bulk)
print_validation_stats("Respired CO2", soilr_ser_resp, npool_ser_resp)

# Feedback model
cat("\nFEEDBACK MODEL:\n")
print_validation_stats("Pool 1 (fast)", soilr_fb_pools[,1], npool_fb_pools[,1])
print_validation_stats("Pool 2 (intermediate)", soilr_fb_pools[,2], npool_fb_pools[,2])
print_validation_stats("Pool 3 (slow)", soilr_fb_pools[,3], npool_fb_pools[,3])
print_validation_stats("Bulk SOM", soilr_fb_bulk, npool_fb_bulk)
print_validation_stats("Respired CO2", soilr_fb_resp, npool_fb_resp)

# =============================================================================
# PART 6: Prepare Data for Visualization
# =============================================================================

# Atmospheric data
atm_df <- data.frame(Year = Fatm[, 1], Delta14C = Fatm[, 2])

# Helper function for pool data (3 pools)
create_pool_comparison_df <- function(years, soilr_pools, npool_pools, model_name) {
  bind_rows(
    data.frame(Year = years, Delta14C = soilr_pools[, 1],
               Pool = "Pool 1 (fast)", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = soilr_pools[, 2],
               Pool = "Pool 2 (intermediate)", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = soilr_pools[, 3],
               Pool = "Pool 3 (slow)", Source = "SoilR", Model = model_name),
    data.frame(Year = years, Delta14C = npool_pools[, 1],
               Pool = "Pool 1 (fast)", Source = "N-Pool", Model = model_name),
    data.frame(Year = years, Delta14C = npool_pools[, 2],
               Pool = "Pool 2 (intermediate)", Source = "N-Pool", Model = model_name),
    data.frame(Year = years, Delta14C = npool_pools[, 3],
               Pool = "Pool 3 (slow)", Source = "N-Pool", Model = model_name)
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
pools_df$Pool <- factor(pools_df$Pool,
                        levels = c("Pool 1 (fast)", "Pool 2 (intermediate)", "Pool 3 (slow)"))

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

# Difference data
diff_df <- bind_rows(
  # Bulk SOM differences
  data.frame(Year = years, Difference = soilr_par_bulk - npool_par_bulk,
             Type = "Bulk SOM", Model = "Parallel"),
  data.frame(Year = years, Difference = soilr_ser_bulk - npool_ser_bulk,
             Type = "Bulk SOM", Model = "Series"),
  data.frame(Year = years, Difference = soilr_fb_bulk - npool_fb_bulk,
             Type = "Bulk SOM", Model = "Feedback"),
  # Respired differences
  data.frame(Year = years, Difference = soilr_par_resp - npool_par_resp,
             Type = "Respired CO2", Model = "Parallel"),
  data.frame(Year = years, Difference = soilr_ser_resp - npool_ser_resp,
             Type = "Respired CO2", Model = "Series"),
  data.frame(Year = years, Difference = soilr_fb_resp - npool_fb_resp,
             Type = "Respired CO2", Model = "Feedback")
)
diff_df$Model <- factor(diff_df$Model, levels = c("Parallel", "Series", "Feedback"))

# =============================================================================
# PART 7: Visualization
# =============================================================================

# Plot 1: Individual pools comparison (one plot per model structure)
p_pools_parallel <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(pools_df, Source == "SoilR", Model == "Parallel"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "black", linewidth = 0.9) +
  geom_line(data = filter(pools_df, Source == "N-Pool", Model == "Parallel"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "red", linewidth = 0.6) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(title = "Parallel", x = NULL, y = expression(Delta^14*C ~ "(‰)"), linetype = NULL) +
  theme_bw() +
  theme(legend.position = "none")

p_pools_series <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(pools_df, Source == "SoilR", Model == "Series"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "black", linewidth = 0.9) +
  geom_line(data = filter(pools_df, Source == "N-Pool", Model == "Series"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "red", linewidth = 0.6) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(title = "Series", x = NULL, y = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "none")

p_pools_feedback <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(pools_df, Source == "SoilR", Model == "Feedback"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "black", linewidth = 0.9) +
  geom_line(data = filter(pools_df, Source == "N-Pool", Model == "Feedback"),
            aes(x = Year, y = Delta14C, linetype = Pool),
            color = "red", linewidth = 0.6) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(title = "Feedback", x = NULL, y = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "right")

# Combined pools plot
p_pools_row <- p_pools_parallel | p_pools_series | p_pools_feedback

# Plot 2: Bulk SOM and Respired comparison
p_bulk <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(bulk_df, Source == "SoilR"),
            aes(x = Year, y = Delta14C, color = Type),
            linewidth = 0.9) +
  geom_line(data = filter(bulk_df, Source == "N-Pool"),
            aes(x = Year, y = Delta14C, color = Type),
            linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~Model, nrow = 1) +
  scale_color_manual(values = c("Bulk SOM" = "blue", "Respired CO2" = "red")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    color = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Plot 3: Difference plots
p_diff <- ggplot(diff_df, aes(x = Year, y = Difference, color = Type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~Model, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Bulk SOM" = "blue", "Respired CO2" = "red")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    x = "Year",
    y = expression(Delta*Delta^14*C ~ "(‰)"),
    color = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Display individual plots
print(p_pools_row)
print(p_bulk)
print(p_diff)

# Combined figure
(p_pools_row) /
  (p_bulk) /
  (p_diff) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Three-Pool Model Validation: SoilR vs N-Pool",
    subtitle = "Black = SoilR (reference) | Red = N-Pool (custom implementation)"
  )

# =============================================================================
# PART 8: Summary Statistics Table
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("SUMMARY: Maximum Absolute Differences (per mil)\n")
cat("=======================================================================\n")
cat("\n")
cat("                    Parallel    Series    Feedback\n")
cat("                    --------    ------    --------\n")

cat(sprintf("Pool 1 (fast)       %8.4f  %8.4f  %8.4f\n",
            max(abs(soilr_par_pools[,1] - npool_par_pools[,1])),
            max(abs(soilr_ser_pools[,1] - npool_ser_pools[,1])),
            max(abs(soilr_fb_pools[,1] - npool_fb_pools[,1]))))

cat(sprintf("Pool 2 (inter.)     %8.4f  %8.4f  %8.4f\n",
            max(abs(soilr_par_pools[,2] - npool_par_pools[,2])),
            max(abs(soilr_ser_pools[,2] - npool_ser_pools[,2])),
            max(abs(soilr_fb_pools[,2] - npool_fb_pools[,2]))))

cat(sprintf("Pool 3 (slow)       %8.4f  %8.4f  %8.4f\n",
            max(abs(soilr_par_pools[,3] - npool_par_pools[,3])),
            max(abs(soilr_ser_pools[,3] - npool_ser_pools[,3])),
            max(abs(soilr_fb_pools[,3] - npool_fb_pools[,3]))))

cat(sprintf("Bulk SOM            %8.4f  %8.4f  %8.4f\n",
            max(abs(soilr_par_bulk - npool_par_bulk)),
            max(abs(soilr_ser_bulk - npool_ser_bulk)),
            max(abs(soilr_fb_bulk - npool_fb_bulk))))

cat(sprintf("Respired CO2        %8.4f  %8.4f  %8.4f\n",
            max(abs(soilr_par_resp - npool_par_resp)),
            max(abs(soilr_ser_resp - npool_ser_resp)),
            max(abs(soilr_fb_resp - npool_fb_resp))))

cat("\n")
cat("Values < 1.0 indicate excellent agreement between implementations.\n")
cat("=======================================================================\n")
