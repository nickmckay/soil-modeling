# =============================================================================
# One-Pool Soil Carbon Model with Radiocarbon (14C) Dynamics
# =============================================================================
#
# This script demonstrates a single-pool soil carbon model with radiocarbon
# tracking using the SoilR package. The one-pool model is the simplest
# representation of soil organic matter dynamics, assuming all soil carbon
# decomposes at a single rate.
#
# Key concepts:
#   - Delta 14C: Radiocarbon expressed as per mil (‰) deviation from a standard
#   - Bomb spike: Nuclear weapons testing (1950s-60s) doubled atmospheric 14C,
#     creating a powerful tracer for carbon cycle studies
#   - Turnover time: Mean time carbon resides in the pool (1/k)
#
# The one-pool model follows:
#   dC/dt = I - k*C
#
# Where:
#   C = carbon stock
#   I = carbon input rate
#   k = decomposition rate constant
#
# References:
#   - Sierra et al. (2017) Models of soil organic matter decomposition
#   - Trumbore (2009) Radiocarbon and Soil Carbon Dynamics
# =============================================================================

library(SoilR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# PART 1: Visualize Atmospheric Radiocarbon Records
# =============================================================================
# SoilR provides two key radiocarbon datasets:
#   - IntCal13: Pre-bomb calibration curve (50,000 years BP to 1950 AD)
#   - Hua2013: Post-bomb atmospheric 14C measurements by hemisphere/zone

# Plot IntCal13 calibration curve with uncertainty bounds (± 1 sigma)
# This shows natural 14C variations over the Holocene and Pleistocene
intcal_df <- data.frame(
  Year = IntCal13$CAL.BP,
  Delta14C = IntCal13$Delta.14C,
  Upper = IntCal13$Delta.14C + IntCal13$Sigma,
  Lower = IntCal13$Delta.14C - IntCal13$Sigma
)

p_intcal <- ggplot(intcal_df, aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3) +
  geom_line(aes(y = Delta14C)) +
  labs(
    title = "IntCal13",
    x = "cal BP",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw()

# Plot Hua et al. (2013) post-bomb curves for different atmospheric zones
# The bomb spike peaked around 1963-1965 and has been declining since
# Different zones show timing differences due to atmospheric mixing
hua_df <- bind_rows(
  data.frame(Year = Hua2013$NHZone1$Year.AD,
             Delta14C = Hua2013$NHZone1$mean.Delta14C,
             Zone = "Northern hemisphere zone 1"),
  data.frame(Year = Hua2013$NHZone2$Year.AD,
             Delta14C = Hua2013$NHZone2$mean.Delta14C,
             Zone = "Northern hemisphere zone 2"),
  data.frame(Year = Hua2013$NHZone3$Year.AD,
             Delta14C = Hua2013$NHZone3$mean.Delta14C,
             Zone = "Northern hemisphere zone 3"),
  data.frame(Year = Hua2013$SHZone12$Year.AD,
             Delta14C = Hua2013$SHZone12$mean.Delta14C,
             Zone = "Southern hemisphere zones 1 and 2"),
  data.frame(Year = Hua2013$SHZone3$Year.AD,
             Delta14C = Hua2013$SHZone3$mean.Delta14C,
             Zone = "Southern hemisphere zone 3")
)

p_hua <- ggplot(hua_df, aes(x = Year, y = Delta14C, color = Zone)) +
  geom_line() +
  labs(
    title = "Hua et al. 2013",
    x = "Year AD",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Display atmospheric radiocarbon records
p_intcal / p_hua

# =============================================================================
# PART 2: Create Combined Atmospheric 14C Record
# =============================================================================
# Bind pre-bomb (IntCal13) and post-bomb (Hua2013) curves into a continuous
# record on the AD time scale. This provides the atmospheric 14C input forcing
# for the soil carbon models.

ad <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone1, time.scale = "AD")

ad_df <- data.frame(Year = ad[, 1], Delta14C = ad[, 2])

p_combined_full <- ggplot(ad_df, aes(x = Year, y = Delta14C)) +
  geom_line() +
  labs(x = "Year", y = expression(Delta^14*C ~ "(‰)")) +
  theme_bw()

p_combined_zoom <- ggplot(ad_df, aes(x = Year, y = Delta14C)) +
  geom_line() +
  geom_vline(xintercept = 1950, linetype = "dashed") +  # 1950: conventional "present" in radiocarbon
  xlim(0, 2010) +
  labs(x = "Year", y = expression(Delta^14*C ~ "(‰)")) +
  theme_bw()

# Display combined atmospheric record
p_combined_full / p_combined_zoom

# =============================================================================
# PART 3: Define Model Parameters
# =============================================================================

# Time vector: simulate from 1901 to 2009 in half-year steps
years <- seq(1901, 2009, by = 0.5)

# Annual litter input (arbitrary units, e.g., g C m^-2 yr^-1)
LitterInput <- 100

# Decomposition rate constants (yr^-1) for different turnover times
# We'll compare pools with different turnover times to see how they respond
# to the bomb spike differently
k_fast <- 1/2      # Turnover time = 2 years (fast cycling)
k_medium <- 1/10   # Turnover time = 10 years (medium cycling)
k_slow <- 1/50     # Turnover time = 50 years (slow cycling)

# Initial carbon stocks (arbitrary units)
# At steady state: C = I/k
C0_fast <- LitterInput / k_fast      # 200 units
C0_medium <- LitterInput / k_medium  # 1000 units
C0_slow <- LitterInput / k_slow      # 5000 units

# Atmospheric radiocarbon forcing (Delta 14C in per mil)
Fatm <- ad[, 1:2]

# =============================================================================
# PART 4: One-Pool Models with Different Turnover Times
# =============================================================================
# In a one-pool model, all soil carbon is treated as a single homogeneous pool.
# The radiocarbon signature of the pool integrates past atmospheric 14C inputs
# weighted by the turnover time.
#
# Structure:  Atmosphere
#                 |
#                 v
#               [Pool]
#                 |
#                 v
#             Respiration
#
# Pools with fast turnover closely track atmospheric 14C.
# Pools with slow turnover show dampened and lagged responses to the bomb spike.

Model_fast <- OnepModel14(
  t = years,
  k = k_fast,
  C0 = C0_fast,
  F0_Delta14C = 0,     # Initial 14C signature (pre-bomb = ~0‰)
  In = LitterInput,
  inputFc = Fatm       # Atmospheric 14C forcing
)

Model_medium <- OnepModel14(
  t = years,
  k = k_medium,
  C0 = C0_medium,
  F0_Delta14C = 0,
  In = LitterInput,
  inputFc = Fatm
)

Model_slow <- OnepModel14(
  t = years,
  k = k_slow,
  C0 = C0_slow,
  F0_Delta14C = 0,
  In = LitterInput,
  inputFc = Fatm
)

# =============================================================================
# PART 5: Extract Radiocarbon Results
# =============================================================================
# getF14R(): 14C signature of respired CO2
# getF14C(): 14C signature of soil carbon stock
# For a one-pool model, these are identical (respiration comes from the same pool)

# Fast turnover model
Fast.R14 <- getF14R(Model_fast)
Fast.C14 <- getF14C(Model_fast)

# Medium turnover model
Medium.R14 <- getF14R(Model_medium)
Medium.C14 <- getF14C(Model_medium)

# Slow turnover model
Slow.R14 <- getF14R(Model_slow)
Slow.C14 <- getF14C(Model_slow)

# =============================================================================
# PART 6: Prepare Data for ggplot2 Visualization
# =============================================================================
# Convert model outputs to tidy data frames for ggplot2

# Atmospheric data for the model period
atm_df <- data.frame(
  Year = Fatm[, 1],
  Delta14C = Fatm[, 2],
  Source = "Atmosphere"
)

# Combine all model results into a single tidy dataframe
pool_df <- bind_rows(
  data.frame(
    Year = years,
    Delta14C = Fast.C14,
    Turnover = "2 years (fast)"
  ),
  data.frame(
    Year = years,
    Delta14C = Medium.C14,
    Turnover = "10 years (medium)"
  ),
  data.frame(
    Year = years,
    Delta14C = Slow.C14,
    Turnover = "50 years (slow)"
  )
)

# Set factor order for legend
pool_df$Turnover <- factor(
  pool_df$Turnover,
  levels = c("2 years (fast)", "10 years (medium)", "50 years (slow)")
)

# =============================================================================
# PART 7: Visualization
# =============================================================================
# Compare how pools with different turnover times respond to the bomb spike
#
# Key observations:
#   - Fast pools (2 yr) closely track atmospheric 14C with minimal lag
#   - Medium pools (10 yr) show a delayed and dampened bomb spike
#   - Slow pools (50 yr) barely register the bomb spike, integrating over decades

# Main comparison plot: All turnover times together
p_comparison <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C, color = "Atmosphere"),
            linewidth = 1) +
  geom_line(data = pool_df, aes(x = Year, y = Delta14C, color = Turnover),
            linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Atmosphere" = "black",
      "2 years (fast)" = "#E41A1C",
      "10 years (medium)" = "#377EB8",
      "50 years (slow)" = "#4DAF4A"
    ),
    breaks = c("Atmosphere", "2 years (fast)", "10 years (medium)", "50 years (slow)")
  ) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "One-Pool Model: Effect of Turnover Time on Radiocarbon Dynamics",
    subtitle = "Faster turnover pools track atmospheric 14C more closely",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    color = "Turnover time"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# Individual faceted plot for detailed comparison
p_faceted <- ggplot() +
  geom_line(data = atm_df, aes(x = Year, y = Delta14C), color = "black") +
  geom_line(data = pool_df, aes(x = Year, y = Delta14C), color = "blue",
            linewidth = 0.8) +
  facet_wrap(~Turnover, ncol = 1) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"))

# Display plots
p_comparison

p_comparison / p_faceted +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(
    title = "One-Pool Radiocarbon Model Comparison",
    subtitle = "Black line = Atmosphere | Colored/Blue lines = Soil pool"
  )
