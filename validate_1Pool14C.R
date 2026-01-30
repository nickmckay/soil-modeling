# =============================================================================
# Validation: One-Pool Model - SoilR vs N-Pool Implementation
# =============================================================================
#
# This script validates the N-pool model implementation against SoilR's
# OnepModel14 by running identical experiments and comparing results visually.
#
# We compare three turnover times (2, 10, and 50 years) using both:
#   - SoilR's OnepModel14 (reference implementation)
#   - Our custom N-pool model (from NPool14C.R)
#
# Success criteria: The two implementations should produce nearly identical
# results (differences < 1 per mil in Delta14C).
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
# Source the N-pool model functions (or define them inline if preferred)

source("NPool14C.R")

# =============================================================================
# PART 2: Set Up Common Parameters
# =============================================================================
# These parameters match exactly those in 1Pool14C.R

# Atmospheric 14C data
ad <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone1, time.scale = "AD")
Fatm <- ad[, 1:2]
atm_data <- prepare_atm_data(Fatm)

# Time vector: simulate from 1901 to 2009 in half-year steps
years <- seq(1901, 2009, by = 0.5)

# Annual litter input
LitterInput <- 100

# Decomposition rate constants for different turnover times
k_fast <- 1/2      # Turnover time = 2 years
k_medium <- 1/10   # Turnover time = 10 years
k_slow <- 1/50     # Turnover time = 50 years

# Initial carbon stocks (at steady state: C = I/k)
C0_fast <- LitterInput / k_fast      # 200 units
C0_medium <- LitterInput / k_medium  # 1000 units
C0_slow <- LitterInput / k_slow      # 5000 units

# Initial 14C (pre-bomb = 0 per mil)
F0_Delta14C <- 0

# =============================================================================
# PART 3: Run SoilR Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Running SoilR One-Pool Models...\n")
cat("=======================================================================\n")

# Fast turnover (2 years)
soilr_fast <- OnepModel14(
  t = years,
  k = k_fast,
  C0 = C0_fast,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  inputFc = Fatm
)

# Medium turnover (10 years)
soilr_medium <- OnepModel14(
  t = years,
  k = k_medium,
  C0 = C0_medium,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  inputFc = Fatm
)

# Slow turnover (50 years)
soilr_slow <- OnepModel14(
  t = years,
  k = k_slow,
  C0 = C0_slow,
  F0_Delta14C = F0_Delta14C,
  In = LitterInput,
  inputFc = Fatm
)

# Extract results
soilr_fast_C14 <- getF14C(soilr_fast)
soilr_medium_C14 <- getF14C(soilr_medium)
soilr_slow_C14 <- getF14C(soilr_slow)

cat("  SoilR models complete.\n")

# =============================================================================
# PART 4: Run N-Pool Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Running N-Pool One-Pool Models...\n")
cat("=======================================================================\n")

# Fast turnover (2 years)
npool_fast <- solve_npool_model(
  times = years,
  k = k_fast,
  C0 = C0_fast,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = 1,
  transfer_matrix = create_parallel_transfers(1),
  F_atm_data = atm_data
)

# Medium turnover (10 years)
npool_medium <- solve_npool_model(
  times = years,
  k = k_medium,
  C0 = C0_medium,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = 1,
  transfer_matrix = create_parallel_transfers(1),
  F_atm_data = atm_data
)

# Slow turnover (50 years)
npool_slow <- solve_npool_model(
  times = years,
  k = k_slow,
  C0 = C0_slow,
  F0 = delta14C_to_F(F0_Delta14C),
  total_input = LitterInput,
  input_fractions = 1,
  transfer_matrix = create_parallel_transfers(1),
  F_atm_data = atm_data
)

# Extract results
npool_fast_C14 <- getF14C_npool(npool_fast)
npool_medium_C14 <- getF14C_npool(npool_medium)
npool_slow_C14 <- getF14C_npool(npool_slow)

cat("  N-Pool models complete.\n")

# =============================================================================
# PART 5: Calculate Differences
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("Validation Results: One-Pool Models\n")
cat("=======================================================================\n")

diff_fast <- soilr_fast_C14 - npool_fast_C14
diff_medium <- soilr_medium_C14 - npool_medium_C14
diff_slow <- soilr_slow_C14 - npool_slow_C14

cat("\n2-Year Turnover (Fast):\n")
cat(sprintf("  Max absolute difference: %.6f ‰\n", max(abs(diff_fast))))
cat(sprintf("  Mean absolute difference: %.6f ‰\n", mean(abs(diff_fast))))
cat(sprintf("  RMSE: %.6f ‰\n", sqrt(mean(diff_fast^2))))

cat("\n10-Year Turnover (Medium):\n")
cat(sprintf("  Max absolute difference: %.6f ‰\n", max(abs(diff_medium))))
cat(sprintf("  Mean absolute difference: %.6f ‰\n", mean(abs(diff_medium))))
cat(sprintf("  RMSE: %.6f ‰\n", sqrt(mean(diff_medium^2))))

cat("\n50-Year Turnover (Slow):\n")
cat(sprintf("  Max absolute difference: %.6f ‰\n", max(abs(diff_slow))))
cat(sprintf("  Mean absolute difference: %.6f ‰\n", mean(abs(diff_slow))))
cat(sprintf("  RMSE: %.6f ‰\n", sqrt(mean(diff_slow^2))))

# =============================================================================
# PART 6: Prepare Data for Visualization
# =============================================================================

# Atmospheric data
atm_df <- data.frame(
  Year = Fatm[, 1],
  Delta14C = Fatm[, 2]
)

# Combine SoilR results
soilr_df <- bind_rows(
  data.frame(Year = years, Delta14C = soilr_fast_C14,
             Turnover = "2 years (fast)", Model = "SoilR"),
  data.frame(Year = years, Delta14C = soilr_medium_C14,
             Turnover = "10 years (medium)", Model = "SoilR"),
  data.frame(Year = years, Delta14C = soilr_slow_C14,
             Turnover = "50 years (slow)", Model = "SoilR")
)

# Combine N-Pool results
npool_df <- bind_rows(
  data.frame(Year = years, Delta14C = npool_fast_C14,
             Turnover = "2 years (fast)", Model = "N-Pool"),
  data.frame(Year = years, Delta14C = npool_medium_C14,
             Turnover = "10 years (medium)", Model = "N-Pool"),
  data.frame(Year = years, Delta14C = npool_slow_C14,
             Turnover = "50 years (slow)", Model = "N-Pool")
)

# Combine all results
all_results_df <- bind_rows(soilr_df, npool_df)
all_results_df$Turnover <- factor(
  all_results_df$Turnover,
  levels = c("2 years (fast)", "10 years (medium)", "50 years (slow)")
)

# Difference data
diff_df <- bind_rows(
  data.frame(Year = years, Difference = diff_fast, Turnover = "2 years (fast)"),
  data.frame(Year = years, Difference = diff_medium, Turnover = "10 years (medium)"),
  data.frame(Year = years, Difference = diff_slow, Turnover = "50 years (slow)")
)
diff_df$Turnover <- factor(
  diff_df$Turnover,
  levels = c("2 years (fast)", "10 years (medium)", "50 years (slow)")
)

# =============================================================================
# PART 7: Visualization
# =============================================================================

# Plot 1: Overlay comparison - SoilR (black) vs N-Pool (red dashed)
p_overlay <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C),
            color = "grey70", linewidth = 0.5) +
  geom_line(data = filter(all_results_df, Model == "SoilR"),
            aes(x = Year, y = Delta14C), color = "black", linewidth = 1) +
  geom_line(data = filter(all_results_df, Model == "N-Pool"),
            aes(x = Year, y = Delta14C), color = "red", linewidth = 0.8,
            linetype = "dashed") +
  facet_wrap(~Turnover, ncol = 1) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "One-Pool Model Validation: SoilR vs N-Pool",
    subtitle = "Black solid = SoilR | Red dashed = N-Pool | Grey = Atmosphere",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"))

# Plot 2: Difference plot
p_diff <- ggplot(diff_df, aes(x = Year, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "darkblue", linewidth = 0.8) +
  facet_wrap(~Turnover, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Difference: SoilR - N-Pool",
    subtitle = "Values near zero indicate good agreement",
    x = "Year",
    y = expression(Delta*Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"))

# Plot 3: Combined comparison with all three turnover times on one plot
p_combined <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C, color = "Atmosphere"),
            linewidth = 0.8) +
  geom_line(data = filter(all_results_df, Model == "SoilR"),
            aes(x = Year, y = Delta14C, color = Turnover),
            linewidth = 1) +
  geom_line(data = filter(all_results_df, Model == "N-Pool"),
            aes(x = Year, y = Delta14C, color = Turnover),
            linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Atmosphere" = "black",
      "2 years (fast)" = "#E41A1C",
      "10 years (medium)" = "#377EB8",
      "50 years (slow)" = "#4DAF4A"
    )
  ) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "One-Pool Model Comparison: All Turnover Times",
    subtitle = "Solid lines = SoilR | Dashed lines = N-Pool",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    color = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Display plots
print(p_overlay)
print(p_diff)
print(p_combined)

# Final combined figure
(p_overlay | p_diff) / p_combined +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "One-Pool Model Validation Summary",
    subtitle = "Comparing SoilR (reference) vs N-Pool (custom implementation)"
  )
