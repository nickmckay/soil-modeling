# Soil Carbon Modeling Project

## Overview

Soil carbon modeling framework in R with radiocarbon (14C) tracking. Models how carbon moves through soil pools and how the atmospheric 14C "bomb spike" (from 1950s-60s nuclear testing) propagates through pools with different turnover times.

## Architecture

### Core solver: `NPool14C.R`
- Generalized N-pool ODE system: `dC/dt = I + A*C` (carbon) and `d(14C)/dt = I*F_atm + A*14C - lambda*14C` (radiocarbon)
- Uses `deSolve::lsoda` for numerical integration
- Key functions: `create_compartmental_matrix()`, `npool_ode()`, `solve_npool_model()`
- Output functions: `getC_npool()`, `getReleaseFlux_npool()`, `getF14_npool()`, `getF14C_npool()`, `getF14R_npool()`
- Utility functions: `delta14C_to_F()`, `F_to_delta14C()`, `prepare_atm_data()`
- Also contains inline validation against SoilR and a 5-pool demo (runs on source)

### Flexible framework: `FlexiblePool14C.R`
- Sources `NPool14C.R` for the core solver
- Builder/config layer for specifying any pool topology via a full n x n flux matrix
- `define_model()` bundles parameters with validation (checks transfer rates don't exceed decomposition rates)
- Convenience constructors: `configure_parallel()`, `configure_series()`, `configure_feedback()`
- `run_model()` wraps `solve_npool_model()` with a config object
- `print_model_summary()` shows readable flux topology
- `plot_model_diagram()` visualizes any model config as a flow diagram (base R graphics, no extra dependencies). Forward transfers curve above pools, backward transfers curve below. Arrow widths scale with rates. Supports custom pool labels and colors.
- `create_flux_matrix(n)` returns an empty n x n matrix for custom topologies

### Earlier standalone models
- `one-pool-model.R` - Initial simple model
- `1Pool14C.R` - One-pool with radiocarbon
- `2Pool14C.R` - Two-pool (parallel, series, feedback)
- `3Pool14C.R` - Three-pool (parallel, series, feedback)

### Validation scripts
- `validate_FlexiblePool14C.R` - Comprehensive: all 9 combinations of {1,2,3} pools x {parallel, series, feedback} vs SoilR. Tests carbon stocks, release fluxes, and Delta14C (bulk SOM, respired CO2, individual pools).
- `validate_1Pool14C.R`, `validate_2Pool14C.R`, `validate_3Pool14C.R` - Earlier per-model-type validations

## Key conventions

- Transfer matrix convention: `T[i,j]` = rate of C flow FROM pool j TO pool i (column = donor, row = receiver)
- Decomposition rates `k` are in yr^-1; turnover time = 1/k
- Delta14C in per mil; fraction modern F = Delta14C/1000 + 1
- Radiocarbon decay constant lambda = log(2)/5730 yr^-1
- Atmospheric 14C data: IntCal13 (pre-bomb) + Hua2013 NHZone1 (post-bomb), bound via `SoilR::bind.C14curves()`

## Dependencies

R packages: `SoilR`, `deSolve`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`

## Running

All scripts are designed to be run from the project root directory. Validation scripts source their dependencies:
```r
source("FlexiblePool14C.R")  # which sources NPool14C.R
```

Run validation: `Rscript validate_FlexiblePool14C.R`

## Known issues

- `NPool14C.R` runs validation and demo code when sourced (side effects). The `FlexiblePool14C.R` framework layer is clean (no side effects beyond sourcing `NPool14C.R`).
- Small numerical differences between our `lsoda` solver and SoilR's solver, especially for fast-cycling pools (2-yr turnover). Bulk SOM < 0.5 per mil, C stocks < 0.002 gC, release fluxes < 0.00002 gC/yr.

## Recent changes

- Fixed bug in `getF14R_npool()`: respiration rate was `k + transfer_out` (wrong) instead of `k - transfer_out`
- Added `getC_npool()` and `getReleaseFlux_npool()` output functions for carbon stocks and release fluxes
- Created `FlexiblePool14C.R` flexible framework with full n x n flux matrix parameterization
- Created `validate_FlexiblePool14C.R` with comprehensive 9-configuration validation including C stocks and fluxes
- Added `plot_model_diagram()` flow visualization function to `FlexiblePool14C.R`
