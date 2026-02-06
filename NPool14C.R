# =============================================================================
# N-Pool Soil Carbon Model with Radiocarbon (14C) Dynamics
# =============================================================================
#
# This script implements a generalized N-pool soil carbon model with radiocarbon
# tracking from first principles. It can represent any number of pools with
# arbitrary transfer structures (parallel, series, feedback, or custom).
#
# The model is built for clarity and explainability rather than computational
# efficiency, using explicit loops and matrix operations.
#
# Mathematical Framework:
# -----------------------
# Carbon dynamics follow a system of linear ODEs:
#
#   dC/dt = I + A * C
#
# Where:
#   C = vector of carbon stocks [C_1, C_2, ..., C_n]
#   I = vector of carbon inputs [I_1, I_2, ..., I_n]
#   A = compartmental matrix (n x n)
#
# The compartmental matrix A encodes:
#   - Diagonal elements A[i,i] = -k_i (decomposition rate of pool i)
#   - Off-diagonal elements A[i,j] = a_ij (transfer rate from pool j to pool i)
#
# Radiocarbon dynamics track 14C content in each pool:
#
#   d(14C_i)/dt = I_i * F_atm + sum_j(A[i,j] * 14C_j) - lambda * 14C_i
#
# Where:
#   14C_i = radiocarbon content in pool i (= C_i * F_i, where F_i is fraction modern)
#   F_atm = atmospheric fraction modern (time-varying)
#   lambda = radiocarbon decay constant (log(2)/5730 ≈ 0.000121 yr^-1)
#
# References:
#   - Sierra et al. (2017) Models of soil organic matter decomposition
#   - Metzler & Sierra (2018) Linear autonomous compartmental models
# =============================================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# For validation against SoilR
library(SoilR)

# =============================================================================
# PART 1: Core Model Functions
# =============================================================================

#' Create a compartmental matrix for an N-pool model
#'
#' @param k Vector of decomposition rates (length n)
#' @param transfer_matrix Matrix of transfer coefficients (n x n), where
#'        transfer_matrix[i,j] is the transfer rate from pool j to pool i.
#'        Diagonal elements are ignored (overwritten by -k).
#' @return Compartmental matrix A (n x n)
create_compartmental_matrix <- function(k, transfer_matrix = NULL) {
  n <- length(k)

  # Initialize matrix
  A <- matrix(0, nrow = n, ncol = n)

  # Add transfer coefficients if provided

if (!is.null(transfer_matrix)) {
    A <- transfer_matrix
  }

  # Set diagonal to negative decomposition rates
  # (this overwrites any diagonal values in transfer_matrix)
  for (i in 1:n) {
    A[i, i] <- -k[i]
  }

  return(A)
}

#' Create transfer matrix for a parallel model (no transfers between pools)
#'
#' @param n Number of pools
#' @return Transfer matrix (n x n) of zeros
create_parallel_transfers <- function(n) {
  return(matrix(0, nrow = n, ncol = n))
}

#' Create transfer matrix for a series model (linear chain)
#'
#' @param k Vector of decomposition rates
#' @param transfer_fractions Vector of fractions transferred to next pool
#'        (length n-1), where transfer_fractions[i] is the fraction of
#'        decomposed C from pool i that goes to pool i+1
#' @return Transfer matrix (n x n)
create_series_transfers <- function(k, transfer_fractions) {
  n <- length(k)
  A <- matrix(0, nrow = n, ncol = n)

  # Transfer from pool i to pool i+1
  for (i in 1:(n - 1)) {
    # a_(i+1, i) = fraction * k_i
    A[i + 1, i] <- transfer_fractions[i] * k[i]
  }

  return(A)
}

#' Create transfer matrix for a feedback model (series with backward transfers)
#'
#' @param k Vector of decomposition rates
#' @param forward_fractions Vector of fractions transferred forward (length n-1)
#' @param backward_fractions Vector of fractions transferred backward (length n-1)
#'        where backward_fractions[i] is the fraction from pool i+1 back to pool i
#' @return Transfer matrix (n x n)
create_feedback_transfers <- function(k, forward_fractions, backward_fractions) {
  n <- length(k)
  A <- matrix(0, nrow = n, ncol = n)

  # Forward transfers: pool i -> pool i+1
  for (i in 1:(n - 1)) {
    A[i + 1, i] <- forward_fractions[i] * k[i]
  }

  # Backward transfers: pool i+1 -> pool i
  for (i in 1:(n - 1)) {
    A[i, i + 1] <- backward_fractions[i] * k[i + 1]
  }

  return(A)
}

#' Create input allocation vector
#'
#' @param total_input Total carbon input rate
#' @param fractions Vector of fractions allocated to each pool (must sum to 1)
#' @return Input vector I
create_input_vector <- function(total_input, fractions) {
  return(total_input * fractions)
}

# =============================================================================
# PART 2: ODE System Definition
# =============================================================================

#' Define the ODE system for carbon and radiocarbon dynamics
#'
#' State vector layout: [C_1, C_2, ..., C_n, 14C_1, 14C_2, ..., 14C_n]
#' First n elements are carbon stocks, next n elements are 14C contents
#'
#' @param t Time
#' @param state State vector (length 2n)
#' @param parms List of parameters
#' @return List containing derivatives
npool_ode <- function(t, state, parms) {
  n <- parms$n
  A <- parms$A
  I <- parms$I
  lambda <- parms$lambda
  F_atm_func <- parms$F_atm_func


  # Extract carbon stocks and 14C contents from state vector
  C <- state[1:n]
  C14 <- state[(n + 1):(2 * n)]

  # Get atmospheric 14C at current time (as fraction modern)
  F_atm <- F_atm_func(t)

  # Carbon dynamics: dC/dt = I + A * C
  dC <- I + A %*% C

  # Radiocarbon dynamics: d(14C)/dt = I * F_atm + A * 14C - lambda * 14C
  # The A matrix applies the same transfers to 14C as to C
  # Additionally, radiocarbon decays with rate lambda
  dC14 <- I * F_atm + A %*% C14 - lambda * C14

  return(list(c(as.vector(dC), as.vector(dC14))))
}

# =============================================================================
# PART 3: Model Solver
# =============================================================================

#' Solve the N-pool model
#'
#' @param times Vector of time points
#' @param k Vector of decomposition rates (length n)
#' @param C0 Initial carbon stocks (length n)
#' @param F0 Initial 14C as fraction modern (length n)
#' @param total_input Total carbon input rate
#' @param input_fractions Fractions of input to each pool (length n, sum to 1)
#' @param transfer_matrix Transfer coefficient matrix (n x n)
#' @param F_atm_data Data frame with columns 'time' and 'F_atm' (fraction modern)
#' @return List with components:
#'         - times: time vector
#'         - C: carbon stocks matrix (times x n)
#'         - C14: 14C content matrix (times x n)
#'         - F: fraction modern matrix (times x n)
#'         - Delta14C: Delta14C matrix (times x n)
solve_npool_model <- function(times, k, C0, F0, total_input, input_fractions,
                               transfer_matrix, F_atm_data) {

  n <- length(k)

  # Build compartmental matrix
  A <- create_compartmental_matrix(k, transfer_matrix)

  # Build input vector
  I <- create_input_vector(total_input, input_fractions)

  # Radiocarbon decay constant (yr^-1)
  # Half-life of 14C is 5730 years
  lambda <- log(2) / 5730

  # Create interpolation function for atmospheric 14C
  # F_atm_data should have fraction modern, not Delta14C
  F_atm_func <- approxfun(
    x = F_atm_data$time,
    y = F_atm_data$F_atm,
    rule = 2  # Use nearest value for extrapolation
  )

  # Initial state: [C_1, ..., C_n, 14C_1, ..., 14C_n]
  # 14C_i = C_i * F_i (14C content = carbon stock * fraction modern)
  C14_0 <- C0 * F0
  state0 <- c(C0, C14_0)

  # Parameters for ODE solver
  parms <- list(
    n = n,
    A = A,
    I = I,
    lambda = lambda,
    F_atm_func = F_atm_func
  )

  # Solve ODE system
  # Using lsoda for stiff/non-stiff automatic switching
  out <- ode(
    y = state0,
    times = times,
    func = npool_ode,
    parms = parms,
    method = "lsoda"
  )

  # Extract results
  C_out <- out[, 2:(n + 1), drop = FALSE]
  C14_out <- out[, (n + 2):(2 * n + 1), drop = FALSE]

  # Calculate fraction modern: F_i = 14C_i / C_i
  F_out <- C14_out / C_out

  # Convert fraction modern to Delta14C (per mil)
  # Delta14C = (F - 1) * 1000
  Delta14C_out <- (F_out - 1) * 1000

  # Store compartmental matrix and input for later calculations
  result <- list(
    times = times,
    C = C_out,
    C14 = C14_out,
    F = F_out,
    Delta14C = Delta14C_out,
    A = A,
    I = I,
    k = k
  )

  return(result)
}

# =============================================================================
# PART 4: Output Functions (analogous to SoilR's getF14, getF14C, getF14R)
# =============================================================================

#' Get Delta14C of each individual pool
#'
#' @param model_output Output from solve_npool_model
#' @return Matrix of Delta14C values (times x n)
getF14_npool <- function(model_output) {
  return(model_output$Delta14C)
}

#' Get Delta14C of bulk soil (stock-weighted mean)
#'
#' @param model_output Output from solve_npool_model
#' @return Vector of bulk soil Delta14C values
getF14C_npool <- function(model_output) {
  C <- model_output$C
  C14 <- model_output$C14

  # Total carbon at each time
  C_total <- rowSums(C)

  # Total 14C at each time
  C14_total <- rowSums(C14)

  # Bulk fraction modern
  F_bulk <- C14_total / C_total

  # Convert to Delta14C
  Delta14C_bulk <- (F_bulk - 1) * 1000

  return(Delta14C_bulk)
}

#' Get Delta14C of respired CO2 (flux-weighted mean)
#'
#' @param model_output Output from solve_npool_model
#' @return Vector of respired CO2 Delta14C values
getF14R_npool <- function(model_output) {
  C <- model_output$C
  C14 <- model_output$C14
  A <- model_output$A
  k <- model_output$k
  n <- length(k)

  # Calculate respiration flux from each pool
  # Respiration = (k_i - sum of transfers out) * C_i
  # But transfers out are encoded in the off-diagonal elements going FROM pool i

  # Sum of transfer coefficients going OUT of each pool
  # These are the column sums of A, excluding the diagonal
  transfer_out <- numeric(n)
  for (j in 1:n) {
    for (i in 1:n) {
      if (i != j) {
        transfer_out[j] <- transfer_out[j] + A[i, j]
      }
    }
  }

  # Respiration rate for each pool = k_i - transfer_out_i
  # (what leaves pool i but doesn't go to another pool)
  resp_rate <- k - transfer_out

  n_times <- nrow(C)
  Delta14C_resp <- numeric(n_times)

  for (t_idx in 1:n_times) {
    # Respiration flux from each pool at this time
    resp_flux <- resp_rate * C[t_idx, ]

    # 14C in respiration from each pool
    # The 14C respired is proportional to the 14C content, not the C content
    resp_flux_14C <- resp_rate * C14[t_idx, ]

    # Total respiration
    total_resp <- sum(resp_flux)
    total_resp_14C <- sum(resp_flux_14C)

    # Flux-weighted fraction modern
    F_resp <- total_resp_14C / total_resp

    # Convert to Delta14C
    Delta14C_resp[t_idx] <- (F_resp - 1) * 1000
  }

  return(Delta14C_resp)
}

#' Get carbon stocks of each individual pool (analogous to SoilR's getC)
#'
#' @param model_output Output from solve_npool_model
#' @return Matrix of carbon stocks (times x n)
getC_npool <- function(model_output) {
  return(model_output$C)
}

#' Get release flux (respiration) from each pool (analogous to SoilR's getReleaseFlux)
#'
#' Respiration from pool j = (k_j - sum of transfer rates out of j) * C_j
#'
#' @param model_output Output from solve_npool_model
#' @return Matrix of release fluxes (times x n)
getReleaseFlux_npool <- function(model_output) {
  C <- model_output$C
  A <- model_output$A
  k <- model_output$k
  n <- length(k)

  # Sum of transfer coefficients going OUT of each pool
  transfer_out <- numeric(n)
  for (j in 1:n) {
    for (i in 1:n) {
      if (i != j) {
        transfer_out[j] <- transfer_out[j] + A[i, j]
      }
    }
  }

  # Respiration rate for each pool
  resp_rate <- k - transfer_out

  # Respiration flux = resp_rate * C (for each time step)
  n_times <- nrow(C)
  release <- matrix(0, nrow = n_times, ncol = n)
  for (j in 1:n) {
    release[, j] <- resp_rate[j] * C[, j]
  }

  return(release)
}

# =============================================================================
# PART 5: Utility Functions
# =============================================================================

#' Convert Delta14C to fraction modern
#'
#' @param Delta14C Delta14C in per mil
#' @return Fraction modern
delta14C_to_F <- function(Delta14C) {
  return(Delta14C / 1000 + 1)
}

#' Convert fraction modern to Delta14C
#'
#' @param F Fraction modern
#' @return Delta14C in per mil
F_to_delta14C <- function(F) {
  return((F - 1) * 1000)
}

#' Prepare atmospheric 14C data for the model
#'
#' @param atm_data Data frame or matrix with time in column 1, Delta14C in column 2
#' @return Data frame with 'time' and 'F_atm' columns
prepare_atm_data <- function(atm_data) {
  data.frame(
    time = atm_data[, 1],
    F_atm = delta14C_to_F(atm_data[, 2])
  )
}

# =============================================================================
# PART 6: Validation Against SoilR Models
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("VALIDATION: Comparing N-Pool Model Against SoilR\n")
cat("=======================================================================\n")

# Load atmospheric 14C data
ad <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone1, time.scale = "AD")
Fatm <- ad[, 1:2]
atm_data <- prepare_atm_data(Fatm)

# Common parameters
years <- seq(1901, 2009, by = 0.5)
LitterInput <- 100

# -----------------------------------------------------------------------------
# Validation 1: One-Pool Model
# -----------------------------------------------------------------------------
cat("\n--- One-Pool Model Validation ---\n")

k1 <- 1/10  # 10-year turnover
C0_1 <- 1000

# SoilR model
soilr_1pool <- OnepModel14(
  t = years,
  k = k1,
  C0 = C0_1,
  F0_Delta14C = 0,
  In = LitterInput,
  inputFc = Fatm
)

# N-pool model (n=1)
npool_1 <- solve_npool_model(
  times = years,
  k = k1,
  C0 = C0_1,
  F0 = delta14C_to_F(0),  # Convert Delta14C to fraction modern
  total_input = LitterInput,
  input_fractions = 1,    # All input to single pool
  transfer_matrix = create_parallel_transfers(1),
  F_atm_data = atm_data
)

# Compare results
soilr_C14 <- getF14C(soilr_1pool)
npool_C14 <- getF14C_npool(npool_1)

max_diff_1 <- max(abs(soilr_C14 - npool_C14))
cat(sprintf("  Max difference in Delta14C: %.4f per mil\n", max_diff_1))
cat(sprintf("  Mean difference: %.4f per mil\n", mean(abs(soilr_C14 - npool_C14))))

# -----------------------------------------------------------------------------
# Validation 2: Two-Pool Parallel Model
# -----------------------------------------------------------------------------
cat("\n--- Two-Pool Parallel Model Validation ---\n")

k2 <- c(1/2, 1/50)
C0_2 <- c(100, 1000)
gam2 <- 0.7  # 70% to fast pool

# SoilR model
soilr_2pool_par <- TwopParallelModel14(
  t = years,
  ks = k2,
  C0 = C0_2,
  F0_Delta14C = c(0, 0),
  In = LitterInput,
  gam = gam2,
  inputFc = Fatm
)

# N-pool model (n=2, parallel)
npool_2_par <- solve_npool_model(
  times = years,
  k = k2,
  C0 = C0_2,
  F0 = c(delta14C_to_F(0), delta14C_to_F(0)),
  total_input = LitterInput,
  input_fractions = c(gam2, 1 - gam2),
  transfer_matrix = create_parallel_transfers(2),
  F_atm_data = atm_data
)

soilr_C14_2par <- getF14C(soilr_2pool_par)
npool_C14_2par <- getF14C_npool(npool_2_par)

max_diff_2par <- max(abs(soilr_C14_2par - npool_C14_2par))
cat(sprintf("  Max difference in bulk Delta14C: %.4f per mil\n", max_diff_2par))

# Also compare individual pools
soilr_pools_2par <- getF14(soilr_2pool_par)
npool_pools_2par <- getF14_npool(npool_2_par)
cat(sprintf("  Max difference in Pool 1: %.4f per mil\n",
            max(abs(soilr_pools_2par[, 1] - npool_pools_2par[, 1]))))
cat(sprintf("  Max difference in Pool 2: %.4f per mil\n",
            max(abs(soilr_pools_2par[, 2] - npool_pools_2par[, 2]))))

# -----------------------------------------------------------------------------
# Validation 3: Two-Pool Series Model
# -----------------------------------------------------------------------------
cat("\n--- Two-Pool Series Model Validation ---\n")

a21_2 <- 0.5 * k2[1]  # 50% transfer from pool 1 to pool 2

# SoilR model
soilr_2pool_ser <- TwopSeriesModel14(
  t = years,
  ks = k2,
  C0 = C0_2,
  F0_Delta14C = c(0, 0),
  In = LitterInput,
  a21 = a21_2,
  inputFc = Fatm
)

# N-pool model (n=2, series)
# Create transfer matrix manually
transfer_2_ser <- matrix(0, 2, 2)
transfer_2_ser[2, 1] <- a21_2  # Transfer from pool 1 to pool 2

npool_2_ser <- solve_npool_model(
  times = years,
  k = k2,
  C0 = C0_2,
  F0 = c(delta14C_to_F(0), delta14C_to_F(0)),
  total_input = LitterInput,
  input_fractions = c(1, 0),  # All input to pool 1 in series model
  transfer_matrix = transfer_2_ser,
  F_atm_data = atm_data
)

soilr_C14_2ser <- getF14C(soilr_2pool_ser)
npool_C14_2ser <- getF14C_npool(npool_2_ser)

max_diff_2ser <- max(abs(soilr_C14_2ser - npool_C14_2ser))
cat(sprintf("  Max difference in bulk Delta14C: %.4f per mil\n", max_diff_2ser))

# -----------------------------------------------------------------------------
# Validation 4: Three-Pool Parallel Model
# -----------------------------------------------------------------------------
cat("\n--- Three-Pool Parallel Model Validation ---\n")

k3 <- c(1/2, 1/10, 1/50)
C0_3 <- c(100, 500, 1000)
gam1_3 <- 0.6
gam2_3 <- 0.2

# SoilR model
soilr_3pool_par <- ThreepParallelModel14(
  t = years,
  ks = k3,
  C0 = C0_3,
  F0_Delta14C = c(0, 0, 0),
  In = LitterInput,
  gam1 = gam1_3,
  gam2 = gam2_3,
  inputFc = Fatm
)

# N-pool model (n=3, parallel)
npool_3_par <- solve_npool_model(
  times = years,
  k = k3,
  C0 = C0_3,
  F0 = rep(delta14C_to_F(0), 3),
  total_input = LitterInput,
  input_fractions = c(gam1_3, gam2_3, 1 - gam1_3 - gam2_3),
  transfer_matrix = create_parallel_transfers(3),
  F_atm_data = atm_data
)

soilr_C14_3par <- getF14C(soilr_3pool_par)
npool_C14_3par <- getF14C_npool(npool_3_par)

max_diff_3par <- max(abs(soilr_C14_3par - npool_C14_3par))
cat(sprintf("  Max difference in bulk Delta14C: %.4f per mil\n", max_diff_3par))

# Compare individual pools
soilr_pools_3par <- getF14(soilr_3pool_par)
npool_pools_3par <- getF14_npool(npool_3_par)
for (i in 1:3) {
  cat(sprintf("  Max difference in Pool %d: %.4f per mil\n", i,
              max(abs(soilr_pools_3par[, i] - npool_pools_3par[, i]))))
}

# -----------------------------------------------------------------------------
# Validation 5: Three-Pool Series Model
# -----------------------------------------------------------------------------
cat("\n--- Three-Pool Series Model Validation ---\n")

a21_3 <- 0.9 * k3[1]
a32_3 <- 0.4 * k3[2]

# SoilR model
soilr_3pool_ser <- ThreepSeriesModel14(
  t = years,
  ks = k3,
  C0 = C0_3,
  F0_Delta14C = c(0, 0, 0),
  In = LitterInput,
  a21 = a21_3,
  a32 = a32_3,
  inputFc = Fatm
)

# N-pool model (n=3, series)
transfer_3_ser <- matrix(0, 3, 3)
transfer_3_ser[2, 1] <- a21_3
transfer_3_ser[3, 2] <- a32_3

npool_3_ser <- solve_npool_model(
  times = years,
  k = k3,
  C0 = C0_3,
  F0 = rep(delta14C_to_F(0), 3),
  total_input = LitterInput,
  input_fractions = c(1, 0, 0),  # All input to pool 1
  transfer_matrix = transfer_3_ser,
  F_atm_data = atm_data
)

soilr_C14_3ser <- getF14C(soilr_3pool_ser)
npool_C14_3ser <- getF14C_npool(npool_3_ser)

max_diff_3ser <- max(abs(soilr_C14_3ser - npool_C14_3ser))
cat(sprintf("  Max difference in bulk Delta14C: %.4f per mil\n", max_diff_3ser))

# Compare respired CO2
soilr_R14_3ser <- getF14R(soilr_3pool_ser)
npool_R14_3ser <- getF14R_npool(npool_3_ser)

max_diff_R_3ser <- max(abs(soilr_R14_3ser - npool_R14_3ser))
cat(sprintf("  Max difference in respired Delta14C: %.4f per mil\n", max_diff_R_3ser))

# -----------------------------------------------------------------------------
# Validation 6: Three-Pool Feedback Model
# -----------------------------------------------------------------------------
cat("\n--- Three-Pool Feedback Model Validation ---\n")

a12_3 <- 0.4 * k3[2]
a23_3 <- 0.7 * k3[3]

# SoilR model
soilr_3pool_fb <- ThreepFeedbackModel14(
  t = years,
  ks = k3,
  C0 = C0_3,
  F0_Delta14C = c(0, 0, 0),
  In = LitterInput,
  a21 = a21_3,
  a12 = a12_3,
  a32 = a32_3,
  a23 = a23_3,
  inputFc = Fatm
)

# N-pool model (n=3, feedback)
transfer_3_fb <- matrix(0, 3, 3)
transfer_3_fb[2, 1] <- a21_3  # Forward: 1 -> 2
transfer_3_fb[1, 2] <- a12_3  # Backward: 2 -> 1
transfer_3_fb[3, 2] <- a32_3  # Forward: 2 -> 3
transfer_3_fb[2, 3] <- a23_3  # Backward: 3 -> 2

npool_3_fb <- solve_npool_model(
  times = years,
  k = k3,
  C0 = C0_3,
  F0 = rep(delta14C_to_F(0), 3),
  total_input = LitterInput,
  input_fractions = c(1, 0, 0),
  transfer_matrix = transfer_3_fb,
  F_atm_data = atm_data
)

soilr_C14_3fb <- getF14C(soilr_3pool_fb)
npool_C14_3fb <- getF14C_npool(npool_3_fb)

max_diff_3fb <- max(abs(soilr_C14_3fb - npool_C14_3fb))
cat(sprintf("  Max difference in bulk Delta14C: %.4f per mil\n", max_diff_3fb))

soilr_R14_3fb <- getF14R(soilr_3pool_fb)
npool_R14_3fb <- getF14R_npool(npool_3_fb)

max_diff_R_3fb <- max(abs(soilr_R14_3fb - npool_R14_3fb))
cat(sprintf("  Max difference in respired Delta14C: %.4f per mil\n", max_diff_R_3fb))

cat("\n=======================================================================\n")
cat("Validation Complete\n")
cat("=======================================================================\n")

# =============================================================================
# PART 7: Visualization of Validation Results
# =============================================================================

# Create comparison plots for the 3-pool series model (representative example)

# Prepare data for plotting
validation_df <- bind_rows(
  # Bulk SOM comparison
  data.frame(
    Year = years,
    Delta14C = soilr_C14_3ser,
    Source = "SoilR",
    Type = "Bulk SOM"
  ),
  data.frame(
    Year = years,
    Delta14C = npool_C14_3ser,
    Source = "N-Pool",
    Type = "Bulk SOM"
  ),
  # Respired CO2 comparison
  data.frame(
    Year = years,
    Delta14C = soilr_R14_3ser,
    Source = "SoilR",
    Type = "Respired CO2"
  ),
  data.frame(
    Year = years,
    Delta14C = npool_R14_3ser,
    Source = "N-Pool",
    Type = "Respired CO2"
  )
)

# Individual pool comparison
pool_validation_df <- bind_rows(
  data.frame(
    Year = rep(years, 3),
    Delta14C = c(soilr_pools_3par[, 1], soilr_pools_3par[, 2], soilr_pools_3par[, 3]),
    Pool = rep(c("Pool 1", "Pool 2", "Pool 3"), each = length(years)),
    Source = "SoilR"
  ),
  data.frame(
    Year = rep(years, 3),
    Delta14C = c(npool_pools_3par[, 1], npool_pools_3par[, 2], npool_pools_3par[, 3]),
    Pool = rep(c("Pool 1", "Pool 2", "Pool 3"), each = length(years)),
    Source = "N-Pool"
  )
)

# Plot 1: Bulk SOM and Respired CO2 comparison
p_validation <- ggplot(validation_df, aes(x = Year, y = Delta14C, color = Source)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Type) +
  scale_color_manual(values = c("SoilR" = "black", "N-Pool" = "red")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Validation: N-Pool Model vs SoilR (3-Pool Series)",
    subtitle = "Red (N-Pool) overlays black (SoilR) - close match indicates validation success",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Plot 2: Individual pool comparison
p_pools_validation <- ggplot(pool_validation_df, aes(x = Year, y = Delta14C, color = Source)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Pool, ncol = 1) +
  scale_color_manual(values = c("SoilR" = "black", "N-Pool" = "red")) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "Validation: Individual Pool Comparison (3-Pool Parallel)",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Display validation plots
p_validation
p_pools_validation

# =============================================================================
# PART 8: Demonstration with Custom N-Pool Configuration
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("DEMONSTRATION: 5-Pool Custom Model\n")
cat("=======================================================================\n")

# Create a 5-pool model with custom structure:
# - Pool 1: Very fast (0.5 yr turnover) - fresh litter
# - Pool 2: Fast (2 yr turnover) - labile SOM
# - Pool 3: Medium (10 yr turnover) - stabilized SOM
# - Pool 4: Slow (50 yr turnover) - mineral-associated
# - Pool 5: Very slow (500 yr turnover) - passive/recalcitrant

k5 <- c(1/0.5, 1/2, 1/10, 1/50, 1/500)
C0_5 <- c(50, 100, 500, 1000, 2000)

# Custom transfer structure: cascading series with some feedback
transfer_5 <- matrix(0, 5, 5)
transfer_5[2, 1] <- 0.8 * k5[1]  # 80% of pool 1 -> pool 2
transfer_5[3, 2] <- 0.6 * k5[2]  # 60% of pool 2 -> pool 3
transfer_5[4, 3] <- 0.4 * k5[3]  # 40% of pool 3 -> pool 4
transfer_5[5, 4] <- 0.2 * k5[4]  # 20% of pool 4 -> pool 5
transfer_5[3, 4] <- 0.1 * k5[4]  # 10% feedback from pool 4 -> pool 3

# All input to pool 1 (litter enters fastest pool)
input_fracs_5 <- c(1, 0, 0, 0, 0)

npool_5 <- solve_npool_model(
  times = years,
  k = k5,
  C0 = C0_5,
  F0 = rep(delta14C_to_F(0), 5),
  total_input = LitterInput,
  input_fractions = input_fracs_5,
  transfer_matrix = transfer_5,
  F_atm_data = atm_data
)

cat("5-Pool model solved successfully.\n")
cat(sprintf("Pool turnover times: %.1f, %.1f, %.1f, %.1f, %.1f years\n",
            1/k5[1], 1/k5[2], 1/k5[3], 1/k5[4], 1/k5[5]))

# Prepare data for plotting
pools_5_df <- data.frame(
  Year = rep(years, 5),
  Delta14C = c(npool_5$Delta14C[, 1], npool_5$Delta14C[, 2],
               npool_5$Delta14C[, 3], npool_5$Delta14C[, 4],
               npool_5$Delta14C[, 5]),
  Pool = rep(c("Pool 1 (0.5 yr)", "Pool 2 (2 yr)", "Pool 3 (10 yr)",
               "Pool 4 (50 yr)", "Pool 5 (500 yr)"), each = length(years))
)
pools_5_df$Pool <- factor(pools_5_df$Pool, levels = unique(pools_5_df$Pool))

bulk_5_df <- data.frame(
  Year = years,
  Delta14C = getF14C_npool(npool_5),
  Source = "Bulk SOM"
)

resp_5_df <- data.frame(
  Year = years,
  Delta14C = getF14R_npool(npool_5),
  Source = "Respired CO2"
)

# Plot individual pools
p_5pool_pools <- ggplot() +
  geom_line(data = data.frame(Year = Fatm[, 1], Delta14C = Fatm[, 2]),
            aes(x = Year, y = Delta14C), color = "black", linewidth = 0.5) +
  geom_line(data = pools_5_df, aes(x = Year, y = Delta14C, color = Pool),
            linewidth = 0.8) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "5-Pool Model: Individual Pool Dynamics",
    subtitle = "Black line = Atmosphere",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Plot bulk and respired
combined_5_df <- bind_rows(bulk_5_df, resp_5_df)

p_5pool_bulk <- ggplot() +
  geom_line(data = data.frame(Year = Fatm[, 1], Delta14C = Fatm[, 2]),
            aes(x = Year, y = Delta14C, color = "Atmosphere"), linewidth = 0.8) +
  geom_line(data = combined_5_df, aes(x = Year, y = Delta14C, color = Source),
            linewidth = 0.8) +
  scale_color_manual(
    values = c("Atmosphere" = "black", "Bulk SOM" = "blue", "Respired CO2" = "red")
  ) +
  coord_cartesian(xlim = c(1940, 2010)) +
  labs(
    title = "5-Pool Model: Bulk SOM vs Respired CO2",
    x = "Year",
    y = expression(Delta^14*C ~ "(‰)")
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Combined display
p_5pool_pools / p_5pool_bulk +
  plot_annotation(
    title = "Custom 5-Pool Radiocarbon Model",
    subtitle = "Demonstrating flexibility of N-Pool framework"
  )
