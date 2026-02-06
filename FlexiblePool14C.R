# =============================================================================
# Flexible N-Pool Soil Carbon Model with Radiocarbon (14C) Dynamics
# =============================================================================
#
# This script provides a flexible framework for specifying soil carbon models
# with any number of pools and arbitrary inter-pool transfer topologies.
#
# The key idea: every model is defined by a full n x n flux matrix where each
# entry specifies the transfer rate between any pair of pools. Specific model
# types (parallel, series, feedback) are just special cases where certain
# fluxes are set to zero.
#
# Model specification:
#   - n pools, each with a decomposition rate k_i
#   - An input allocation vector (fraction of litter input to each pool)
#   - A transfer matrix T (n x n) where T[i,j] = rate of C flow from pool j
#     to pool i. Diagonal elements are ignored (overwritten by -k_i).
#
# This sources the core solver from NPool14C.R and adds a builder/config layer.
#
# =============================================================================

# Load core model functions (ODE solver, compartmental matrix, output functions)
source("NPool14C.R")

# =============================================================================
# Model Configuration Functions
# =============================================================================

#' Create an empty flux matrix for an N-pool model
#'
#' Returns an n x n matrix of zeros that can be filled in to specify
#' transfer rates between any pair of pools.
#'
#' Convention: flux_matrix[i, j] = rate of carbon transfer from pool j to pool i
#'
#' @param n Number of pools
#' @return n x n matrix of zeros
create_flux_matrix <- function(n) {
  m <- matrix(0, nrow = n, ncol = n)
  rownames(m) <- paste0("to_pool_", 1:n)
  colnames(m) <- paste0("from_pool_", 1:n)
  return(m)
}

#' Define a model configuration with validation
#'
#' Bundles all model parameters into a named list and validates that the
#' transfer rates are physically consistent (no pool transfers out more
#' carbon than it decomposes).
#'
#' @param k Vector of decomposition rates (length n), or turnover times if
#'        as_turnover_times = TRUE
#' @param input_fractions Vector of fractions of litter input to each pool
#'        (length n, must sum to 1)
#' @param flux_matrix n x n transfer matrix (T[i,j] = rate from pool j to pool i)
#' @param as_turnover_times If TRUE, interpret k as turnover times (years) and
#'        convert to rates internally
#' @return Named list (model config) with components: n, k, input_fractions,
#'         flux_matrix, turnover_times
define_model <- function(k, input_fractions, flux_matrix,
                          as_turnover_times = FALSE) {

  # Convert turnover times to rates if needed
  if (as_turnover_times) {
    turnover_times <- k
    k <- 1 / k
  } else {
    turnover_times <- 1 / k
  }

  n <- length(k)

  # --- Validation ---

  # Check dimensions
  if (length(input_fractions) != n) {
    stop(sprintf("input_fractions length (%d) must match number of pools (%d)",
                 length(input_fractions), n))
  }
  if (!all(dim(flux_matrix) == c(n, n))) {
    stop(sprintf("flux_matrix must be %d x %d", n, n))
  }

  # Check input fractions sum to 1
  if (abs(sum(input_fractions) - 1) > 1e-10) {
    stop(sprintf("input_fractions must sum to 1 (got %.6f)", sum(input_fractions)))
  }

  # Check that transfer rates out of each pool don't exceed decomposition rate
  # Column sum of off-diagonal elements = total transfer rate out of pool j
  for (j in 1:n) {
    total_transfer_out <- sum(flux_matrix[, j]) - flux_matrix[j, j]
    if (total_transfer_out > k[j] + 1e-10) {
      warning(sprintf(
        "Pool %d: total transfer out (%.4f) exceeds decomposition rate (%.4f). ",
        j, total_transfer_out, k[j]),
        "This means the pool transfers more than it decomposes, which is ",
        "physically inconsistent.")
    }
  }

  # Zero out diagonal of flux matrix (these are set by k, not by user)
  diag(flux_matrix) <- 0

  config <- list(
    n = n,
    k = k,
    turnover_times = turnover_times,
    input_fractions = input_fractions,
    flux_matrix = flux_matrix
  )

  class(config) <- "soil_model_config"
  return(config)
}

#' Configure a parallel model (no inter-pool transfers)
#'
#' All pools receive direct atmospheric input; no carbon is transferred
#' between pools.
#'
#' @param k Vector of decomposition rates (length n)
#' @param input_fractions Vector of fractions of litter input to each pool
#' @param as_turnover_times If TRUE, interpret k as turnover times
#' @return Model config
configure_parallel <- function(k, input_fractions, as_turnover_times = FALSE) {
  n <- length(k)
  flux_matrix <- create_flux_matrix(n)
  define_model(k, input_fractions, flux_matrix,
               as_turnover_times = as_turnover_times)
}

#' Configure a series model (linear chain)
#'
#' All input enters pool 1. Carbon flows sequentially: pool 1 -> pool 2 -> ...
#' transfer_fractions[i] is the fraction of decomposed C from pool i that
#' passes to pool i+1 (the rest is respired).
#'
#' @param k Vector of decomposition rates (length n)
#' @param transfer_fractions Vector of fractions transferred forward (length n-1)
#' @param as_turnover_times If TRUE, interpret k as turnover times
#' @return Model config
configure_series <- function(k, transfer_fractions, as_turnover_times = FALSE) {
  if (as_turnover_times) {
    k_rates <- 1 / k
  } else {
    k_rates <- k
  }
  n <- length(k_rates)

  if (n == 1) {
    # Single pool: no transfers possible
    flux_matrix <- create_flux_matrix(1)
    return(define_model(k, input_fractions = 1, flux_matrix,
                        as_turnover_times = as_turnover_times))
  }

  if (length(transfer_fractions) != n - 1) {
    stop(sprintf("transfer_fractions length (%d) must be n-1 (%d)",
                 length(transfer_fractions), n - 1))
  }

  flux_matrix <- create_flux_matrix(n)
  for (i in 1:(n - 1)) {
    flux_matrix[i + 1, i] <- transfer_fractions[i] * k_rates[i]
  }

  input_fractions <- c(1, rep(0, n - 1))

  define_model(k, input_fractions, flux_matrix,
               as_turnover_times = as_turnover_times)
}

#' Configure a feedback model (series with backward transfers)
#'
#' All input enters pool 1. Carbon flows forward (pool i -> pool i+1) and
#' backward (pool i+1 -> pool i) between adjacent pools.
#'
#' @param k Vector of decomposition rates (length n)
#' @param forward_fractions Vector of fractions transferred forward (length n-1)
#' @param backward_fractions Vector of fractions transferred backward (length n-1)
#'        where backward_fractions[i] is the fraction from pool i+1 back to pool i
#' @param as_turnover_times If TRUE, interpret k as turnover times
#' @return Model config
configure_feedback <- function(k, forward_fractions, backward_fractions,
                                as_turnover_times = FALSE) {
  if (as_turnover_times) {
    k_rates <- 1 / k
  } else {
    k_rates <- k
  }
  n <- length(k_rates)

  if (n == 1) {
    flux_matrix <- create_flux_matrix(1)
    return(define_model(k, input_fractions = 1, flux_matrix,
                        as_turnover_times = as_turnover_times))
  }

  if (length(forward_fractions) != n - 1) {
    stop(sprintf("forward_fractions length (%d) must be n-1 (%d)",
                 length(forward_fractions), n - 1))
  }
  if (length(backward_fractions) != n - 1) {
    stop(sprintf("backward_fractions length (%d) must be n-1 (%d)",
                 length(backward_fractions), n - 1))
  }

  flux_matrix <- create_flux_matrix(n)

  # Forward transfers: pool i -> pool i+1
  for (i in 1:(n - 1)) {
    flux_matrix[i + 1, i] <- forward_fractions[i] * k_rates[i]
  }

  # Backward transfers: pool i+1 -> pool i
  for (i in 1:(n - 1)) {
    flux_matrix[i, i + 1] <- backward_fractions[i] * k_rates[i + 1]
  }

  input_fractions <- c(1, rep(0, n - 1))

  define_model(k, input_fractions, flux_matrix,
               as_turnover_times = as_turnover_times)
}

#' Run a model from a configuration
#'
#' @param config Model config from define_model or configure_* functions
#' @param times Vector of time points
#' @param C0 Initial carbon stocks (length n)
#' @param F0 Initial 14C as fraction modern (length n)
#' @param total_input Total carbon input rate
#' @param F_atm_data Data frame with 'time' and 'F_atm' columns
#' @return Output from solve_npool_model
run_model <- function(config, times, C0, F0, total_input, F_atm_data) {
  solve_npool_model(
    times = times,
    k = config$k,
    C0 = C0,
    F0 = F0,
    total_input = total_input,
    input_fractions = config$input_fractions,
    transfer_matrix = config$flux_matrix,
    F_atm_data = F_atm_data
  )
}

#' Visualize the model configuration as a flow diagram
#'
#' Draws pools as boxes with arrows showing atmospheric inputs, inter-pool
#' transfers, and respiration outputs. Arrow widths are proportional to rates.
#' Works for any number of pools and any transfer topology.
#'
#' Forward transfers (left to right) curve above the pools; backward transfers
#' (right to left) curve below. Non-adjacent transfers arc higher to clear
#' intermediate pools.
#'
#' @param config Model config from define_model or configure_* functions
#' @param total_input Optional total litter input rate for labeling (default 100)
#' @param title Optional plot title
#' @param pool_labels Optional vector of pool names (default "Pool 1", "Pool 2", ...)
#' @param col_pools Color for pool boxes (default "cornflowerblue")
#' @param col_input Color for input arrows (default "forestgreen")
#' @param col_transfer Color for transfer arrows (default "darkorange")
#' @param col_resp Color for respiration arrows (default "firebrick")
plot_model_diagram <- function(config, total_input = 100, title = NULL,
                                pool_labels = NULL,
                                col_pools = "cornflowerblue",
                                col_input = "forestgreen",
                                col_transfer = "darkorange",
                                col_resp = "firebrick") {
  n <- config$n
  k <- config$k
  flux <- config$flux_matrix
  input_frac <- config$input_fractions
  tau <- config$turnover_times

  if (is.null(pool_labels)) {
    pool_labels <- paste0("Pool ", 1:n)
  }
  if (is.null(title)) {
    title <- sprintf("Soil Carbon Model: %d pool(s)", n)
  }

  # --- Layout parameters ---
  box_h <- 0.06                    # box half-height
  pool_y <- 0.48                   # vertical center of pool row
  atm_y <- 0.92                    # atmosphere label
  resp_y <- 0.08                   # respiration label

  # Box width adapts to number of pools and longest label
  max_label_nchar <- max(nchar(pool_labels))
  box_w <- min(0.12, 0.8 / (n + 0.5))
  box_w <- max(box_w, max_label_nchar * 0.008 + 0.02)  # ensure labels fit

  # Spread pools evenly
  if (n == 1) {
    pool_x <- 0.5
  } else {
    margin <- box_w + 0.06
    pool_x <- seq(margin, 1 - margin, length.out = n)
  }

  # Gap between adjacent pool centers (used for arrow curve scaling)
  pool_gap <- if (n > 1) pool_x[2] - pool_x[1] else 0.3

  # Arrow line width scaling: map fraction of k (0 to 1) to lwd [1, 5]
  scale_lwd <- function(frac) {
    pmax(1, 1 + 4 * pmin(frac, 1))
  }

  # --- Set up plot ---
  old_par <- par(mar = c(1, 1, 3, 1), xpd = TRUE)
  on.exit(par(old_par))

  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = NA,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = title)

  # --- Helper: draw a curved arrow between two pool positions ---
  draw_curved_arrow <- function(x_from, y_from, x_to, y_to,
                                 curve_height, lwd, col, label = NULL) {
    n_pts <- 60
    t_seq <- seq(0, 1, length.out = n_pts)

    # Control point for quadratic bezier
    mid_x <- (x_from + x_to) / 2
    mid_y <- (y_from + y_to) / 2 + curve_height

    # Quadratic bezier
    bx <- (1 - t_seq)^2 * x_from + 2 * (1 - t_seq) * t_seq * mid_x +
          t_seq^2 * x_to
    by <- (1 - t_seq)^2 * y_from + 2 * (1 - t_seq) * t_seq * mid_y +
          t_seq^2 * y_to

    lines(bx, by, lwd = lwd, col = col)

    # Arrowhead at the end
    arrows(bx[n_pts - 2], by[n_pts - 2], bx[n_pts], by[n_pts],
           length = 0.08, lwd = lwd, col = col)

    # Label at the peak of the curve
    if (!is.null(label)) {
      label_offset <- ifelse(curve_height > 0, 0.025, -0.025)
      text(mid_x, mid_y + label_offset, label, cex = 0.7, col = col)
    }
  }

  # --- Draw pool boxes ---
  for (i in 1:n) {
    rect(pool_x[i] - box_w, pool_y - box_h,
         pool_x[i] + box_w, pool_y + box_h,
         col = adjustcolor(col_pools, alpha.f = 0.15),
         border = col_pools, lwd = 2)
    text(pool_x[i], pool_y + 0.015,
         pool_labels[i], cex = 0.85, font = 2)
    text(pool_x[i], pool_y - 0.022,
         bquote(tau ~ "=" ~ .(format(tau[i], digits = 3)) ~ "yr"),
         cex = 0.65, col = "grey30")
  }

  # --- Draw atmosphere label ---
  text(0.5, atm_y, "Atmosphere", cex = 1.1, font = 2, col = col_input)

  # --- Draw input arrows (atmosphere -> pools) ---
  for (i in 1:n) {
    if (input_frac[i] > 0) {
      lwd <- scale_lwd(input_frac[i])

      arrows(pool_x[i], atm_y - 0.05,
             pool_x[i], pool_y + box_h + 0.01,
             length = 0.10, lwd = lwd, col = col_input)

      # Label to the right of the arrow
      label_y <- (atm_y - 0.05 + pool_y + box_h) / 2
      text(pool_x[i] + box_w * 0.7, label_y,
           sprintf("%.0f%%", input_frac[i] * 100),
           cex = 0.7, col = col_input, adj = 0)
    }
  }

  # --- Draw respiration arrows (pools -> out) ---
  text(0.5, resp_y - 0.05, "Respiration", cex = 1.0, font = 2, col = col_resp)

  for (j in 1:n) {
    transfer_out_j <- sum(flux[, j]) - flux[j, j]
    frac_resp <- 1 - transfer_out_j / k[j]

    if (frac_resp > 0.001) {
      lwd <- scale_lwd(frac_resp)
      arrows(pool_x[j], pool_y - box_h - 0.01,
             pool_x[j], resp_y + 0.01,
             length = 0.10, lwd = lwd, col = col_resp)

      label_y <- (pool_y - box_h + resp_y) / 2
      text(pool_x[j] + box_w * 0.7, label_y,
           sprintf("%.0f%%", frac_resp * 100),
           cex = 0.7, col = col_resp, adj = 0)
    }
  }

  # --- Draw inter-pool transfer arrows ---
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && flux[i, j] > 0) {
        frac_of_k <- flux[i, j] / k[j]
        lwd <- scale_lwd(frac_of_k)

        # Number of pools spanned (for scaling arc height)
        span <- abs(i - j)

        # Forward (left to right): depart from top-right of donor, arrive
        # at top-left of receiver, curve above
        # Backward (right to left): depart from bottom-left of donor,
        # arrive at bottom-right of receiver, curve below
        if (j < i) {
          # Forward transfer (j -> i, left to right)
          x0 <- pool_x[j] + box_w
          y0 <- pool_y + box_h * 0.5
          x1 <- pool_x[i] - box_w
          y1 <- pool_y + box_h * 0.5
          curve_h <- 0.08 + 0.05 * (span - 1)
        } else {
          # Backward transfer (j -> i, right to left)
          x0 <- pool_x[j] - box_w
          y0 <- pool_y - box_h * 0.5
          x1 <- pool_x[i] + box_w
          y1 <- pool_y - box_h * 0.5
          curve_h <- -(0.08 + 0.05 * (span - 1))
        }

        label <- sprintf("%.0f%%", frac_of_k * 100)
        draw_curved_arrow(x0, y0, x1, y1, curve_h, lwd,
                           col_transfer, label)
      }
    }
  }

  # --- Legend ---
  legend(0.0, 0.06,
         legend = c("Atmospheric input", "Inter-pool transfer", "Respiration"),
         col = c(col_input, col_transfer, col_resp),
         lwd = 2, cex = 0.7, bty = "n", horiz = TRUE)
}

#' Print a readable summary of a model configuration
#'
#' @param config Model config
print_model_summary <- function(config) {
  n <- config$n

  cat(sprintf("Soil Carbon Model: %d pool(s)\n", n))
  cat(sprintf("  %-10s %-15s %-15s %-15s\n",
              "Pool", "Turnover (yr)", "k (yr^-1)", "Input fraction"))
  cat(sprintf("  %-10s %-15s %-15s %-15s\n",
              "----", "-------------", "---------", "--------------"))
  for (i in 1:n) {
    cat(sprintf("  %-10d %-15.2f %-15.4f %-15.3f\n",
                i, config$turnover_times[i], config$k[i],
                config$input_fractions[i]))
  }

  # Print non-zero transfers
  has_transfers <- FALSE
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && config$flux_matrix[i, j] > 0) {
        if (!has_transfers) {
          cat("\n  Inter-pool transfers:\n")
          has_transfers <- TRUE
        }
        # Express as fraction of donor pool's k
        frac_of_k <- config$flux_matrix[i, j] / config$k[j]
        cat(sprintf("    Pool %d -> Pool %d: %.4f yr^-1 (%.0f%% of k_%d)\n",
                    j, i, config$flux_matrix[i, j], frac_of_k * 100, j))
      }
    }
  }
  if (!has_transfers) {
    cat("\n  No inter-pool transfers (parallel structure)\n")
  }

  # Calculate and display respiration fractions
  cat("\n  Carbon fate per pool (fraction of decomposed C):\n")
  for (j in 1:n) {
    total_transfer_out <- sum(config$flux_matrix[, j]) - config$flux_matrix[j, j]
    frac_transferred <- total_transfer_out / config$k[j]
    frac_respired <- 1 - frac_transferred
    cat(sprintf("    Pool %d: %.0f%% respired, %.0f%% transferred\n",
                j, frac_respired * 100, frac_transferred * 100))
  }
}
