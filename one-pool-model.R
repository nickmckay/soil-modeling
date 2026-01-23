# One-Pool Soil Carbon Model
#
# This script implements a simple one-pool soil carbon model that simulates
# the dynamics of soil organic carbon over time. The model assumes carbon
# enters the soil at a constant rate and decomposes following first-order
# kinetics (decomposition is proportional to the current carbon stock).
#
# The governing differential equation is:
#   dC/dt = In - k * C
# where:
#   C  = carbon stock (mass)
#   In = carbon input rate (mass/time)
#   k  = decomposition rate constant (1/time)
#
# At steady state (dC/dt = 0), the equilibrium carbon stock = In/k

library(SoilR)

# Define time parameters (in years)
t_start = 0
t_end = 10
t_step = 1/12
t = seq(t_start, t_end, t_step)

# Redefine time vector with coarser step for comparison
t = seq(0, 10, by = 0.2)

# Model parameters
k = 0.8   # Decomposition rate constant (1/year) - higher = faster decay
C0 = 100  # Initial carbon stock (kg)
In = 30   # Carbon input rate (kg/year) - e.g., from litter, root inputs

# -----------------------------------------------------------------------------
# SoilR Package Implementation
# -----------------------------------------------------------------------------
# OnepModel() creates a one-pool model object using the SoilR package
# getC() extracts the carbon stock time series from the model
Model1 = OnepModel(t, k, C0, In)
Ct1 = getC(Model1)

# -----------------------------------------------------------------------------
# Custom Implementation (for comparison/validation)
# -----------------------------------------------------------------------------
# Manual implementation using a simple forward method to solve the ODE
# This helps verify understanding of the model and validate against SoilR
cstock <- c()
for (ti in seq_along(t)) {
  if (ti == 1) {
    # Initialize with starting carbon stock
    cstock[ti] <- C0
  } else {
    # Calculate time step size
    timestep <- t[ti] - t[ti - 1]
    # Calculate change in carbon: input minus decomposition losses
    # Note: This uses a simplified discrete approximation
    dCdT <- In - cstock[ti - 1] * (k * timestep)
    # Update carbon stock
    cstock[ti] <- cstock[ti - 1] + dCdT
  }
}

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
# Plot both model outputs to compare SoilR (black) vs custom implementation (red)
# Deviations indicate differences in numerical methods or implementation errors
plot(t, Ct1, type = "l", ylab = "Carbon stocks (kg)", xlab = "Time (years)")
lines(t, cstock, col = "red")
