# =============================================================================
# Two-Pool Soil Carbon Model with Radiocarbon (14C) Dynamics
# =============================================================================
#
# This script demonstrates two-pool soil carbon models with radiocarbon
# tracking using the SoilR package. It compares how different model structures
# (parallel, series, feedback) affect the radiocarbon signature of soil organic
# matter (SOM) and respired CO2.
#
# Key concepts:
#   - Delta 14C: Radiocarbon expressed as per mil (‰) deviation from a standard
#   - Bomb spike: Nuclear weapons testing (1950s-60s) doubled atmospheric 14C,
#     creating a powerful tracer for carbon cycle studies
#   - Transit time: Time carbon spends from entry to exit from the system
#   - System age: Age of carbon currently residing in the system
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

# Decomposition rate constants (yr^-1) for two pools
# k = 1/turnover_time, so:
#   Pool 1 (fast): k1 = 1/2 yr^-1 → turnover time = 2 years
#   Pool 2 (slow): k2 = 1/50 yr^-1 → turnover time = 50 years
k <- c(k1 = 1/2, k2 = 1/50)

# Initial carbon stocks for each pool (arbitrary units)
C0 <- c(100, 1000)

# Atmospheric radiocarbon forcing (Delta 14C in per mil)
Fatm <- ad[, 1:2]

# =============================================================================
# PART 4: Two-Pool Parallel Model
# =============================================================================
# In a parallel model, each pool receives inputs directly from the atmosphere
# (via litter) and decomposes independently. There is no transfer between pools.
#
# Structure:    Atmosphere
#                 /    \
#                v      v
#              [P1]    [P2]  (independent pools)
#                \      /
#                 v    v
#              Respiration
#
# gamma: Fraction of litter input allocated to pool 1
# (Pool 2 receives 1 - gamma)

Parallel <- TwopParallelModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = c(0, 0),  # Initial 14C signature of each pool (pre-bomb = ~0‰)
  In = LitterInput,
  gam = 0.7,              # 70% of input goes to fast pool
  inputFc = Fatm          # Atmospheric 14C forcing
)

# =============================================================================
# PART 5: Two-Pool Series Model
# =============================================================================
# In a series model, carbon flows sequentially through the pools:
# fast → slow. This represents progressive stabilization of organic matter
# as it ages.
#
# Structure:  Atmosphere
#                 |
#                 v
#               [P1] → [P2]
#                 |      |
#                 v      v
#               Respiration (from each pool)
#
# a21: Transfer coefficient from pool 1 to pool 2 (fraction of k1 * C1)

a21 <- 0.5 * k[1]  # 50% of decomposed fast pool C transfers to slow pool

Series <- TwopSeriesModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = c(0, 0),
  In = LitterInput,
  a21 = a21,
  inputFc = Fatm
)

# =============================================================================
# PART 6: Two-Pool Feedback Model
# =============================================================================
# The feedback model extends the series model by allowing carbon to cycle
# back from the slow to fast pool. This represents processes like microbial
# recycling or physical destabilization of protected carbon.
#
# Structure:  Atmosphere
#                 |
#                 v
#               [P1] ⇄ [P2]
#                 |      |
#                 v      v
#               Respiration (from each pool)
#
# a12: Backward transfer coefficient (feedback from slow to fast)

a12 <- 0.3 * k[2]  # 30% of decomposed slow pool C returns to fast pool

Feedback <- TwopFeedbackModel14(
  t = years,
  ks = k,
  C0 = C0,
  F0_Delta14C = c(0, 0),
  In = LitterInput,
  a21 = a21,
  a12 = a12,
  inputFc = Fatm
)

# =============================================================================
# PART 7: Extract Radiocarbon Results
# =============================================================================
# getF14R(): 14C signature of respired CO2 (flux-weighted mean)
# getF14C(): 14C signature of bulk SOM (stock-weighted mean of all pools)
# getF14():  14C signature of each individual pool

# Parallel model results
P.R14m <- getF14R(Parallel)  # Respired 14C
P.C14m <- getF14C(Parallel)  # Bulk SOM 14C
P.C14t <- getF14(Parallel)   # Individual pool 14C

# Series model results
S.R14m <- getF14R(Series)
S.C14m <- getF14C(Series)
S.C14t <- getF14(Series)

# Feedback model results
F.R14m <- getF14R(Feedback)
F.C14m <- getF14C(Feedback)
F.C14t <- getF14(Feedback)

# =============================================================================
# PART 8: Prepare Data for ggplot2 Visualization
# =============================================================================
# Convert model outputs to tidy data frames for ggplot2

# Atmospheric data for the model period
atm_df <- data.frame(
  Year = Fatm[, 1],
  Delta14C = Fatm[, 2],
  Source = "Atmosphere"
)

# Helper function to create tidy dataframe for individual pool results
create_pool_df <- function(years, pool_data, model_name) {
  data.frame(
    Year = rep(years, 2),
    Delta14C = c(pool_data[, 1], pool_data[, 2]),
    Pool = rep(c("Pool 1 (fast)", "Pool 2 (slow)"), each = length(years)),
    Model = model_name
  )
}

# Helper function to create tidy dataframe for bulk/respired results
create_bulk_df <- function(years, bulk_data, resp_data, model_name) {
  data.frame(
    Year = rep(years, 2),
    Delta14C = c(bulk_data, resp_data),
    Source = rep(c("Bulk SOM", "Respired C"), each = length(years)),
    Model = model_name
  )
}

# Combine individual pool data from all models
pools_df <- bind_rows(
  create_pool_df(years, P.C14t, "Parallel"),
  create_pool_df(years, S.C14t, "Series"),
  create_pool_df(years, F.C14t, "Feedback")
)
pools_df$Model <- factor(pools_df$Model, levels = c("Parallel", "Series", "Feedback"))

# Combine bulk SOM and respired C data from all models
bulk_df <- bind_rows(
  create_bulk_df(years, P.C14m, P.R14m, "Parallel"),
  create_bulk_df(years, S.C14m, S.R14m, "Series"),
  create_bulk_df(years, F.C14m, F.R14m, "Feedback")
)
bulk_df$Model <- factor(bulk_df$Model, levels = c("Parallel", "Series", "Feedback"))

# =============================================================================
# PART 9: Visualization
# =============================================================================
# Create a 3x2 panel comparing all three model structures
# Left column: Individual pool 14C dynamics
# Right column: Bulk SOM vs respired 14C
#
# Key observations to look for:
#   - Fast pools closely track atmospheric 14C
#   - Slow pools show delayed and dampened bomb spike response
#   - Respired 14C (red) reflects the flux-weighted age of decomposing C
#   - Bulk SOM 14C (blue) reflects the stock-weighted mean age

# Left column: Individual pool trajectories for each model
p_pools <- ggplot() +
  # Atmosphere as background reference
  geom_line(data = atm_df, aes(x = Year, y = Delta14C), color = "black") +
  # Individual pools
  geom_line(data = pools_df, aes(x = Year, y = Delta14C, linetype = Pool), color = "blue") +
  facet_wrap(~Model, ncol = 1) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)"),
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Right column: Bulk SOM vs respired CO2 for each model
p_bulk <- ggplot() +
  # Atmosphere as background reference
  geom_line(data = atm_df, aes(x = Year, y = Delta14C, color = "Atmosphere")) +
  # Bulk SOM and respired C
  geom_line(data = bulk_df, aes(x = Year, y = Delta14C, color = Source)) +
  facet_wrap(~Model, ncol = 1) +
  scale_color_manual(
    values = c("Atmosphere" = "black", "Bulk SOM" = "blue", "Respired C" = "red"),
    breaks = c("Atmosphere", "Bulk SOM", "Respired C")
  ) +
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

# Combine into final 2-column layout using patchwork
p_pools + p_bulk +
  plot_annotation(
    title = "Two-Pool Radiocarbon Model Comparison",
    subtitle = "Left: Individual pool dynamics | Right: Bulk SOM vs Respired CO2"
  )
