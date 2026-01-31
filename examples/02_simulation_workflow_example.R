# Example 2: Seroincidence Simulation Workflow
# This script demonstrates how to run seroincidence simulations
# using the shigella package functions

library(shigella)
library(tibble)
library(dplyr)
library(serodynamics)

# -----------------------------------------------------------------------------
# 1. Create Mock Antibody Decay Curve Parameters
# -----------------------------------------------------------------------------
# Note: Real curve parameters will be obtained from MCMC fitting of
# longitudinal data once real data is available

set.seed(456)

# Generate mock MCMC samples of curve parameters
# These represent posterior distributions from Bayesian curve fitting
n_mcmc <- 100

mock_curve_params <- tibble(
  antigen_iso = rep("IgG", n_mcmc),
  iter = 1:n_mcmc,
  chain = rep(1, n_mcmc),
  # Antibody peak parameters
  y0 = rnorm(n_mcmc, 2, 0.3),          # Baseline log(antibody)
  y1 = rnorm(n_mcmc, 5, 0.5),          # Peak log(antibody)
  t1 = rnorm(n_mcmc, 30, 5),           # Days to peak
  # Decay parameters
  alpha = rnorm(n_mcmc, 3, 0.5),       # Short-term decay rate
  r = abs(rnorm(n_mcmc, 0.01, 0.003))  # Long-term decay rate
)

# Set class and attributes for serodynamics
class(mock_curve_params) <- c("curve_params", class(mock_curve_params))
attr(mock_curve_params, "antigen_isos") <- c("IgG")

cat("Mock curve parameters created:\n")
print(head(mock_curve_params))

# -----------------------------------------------------------------------------
# 2. Define Noise Parameters
# -----------------------------------------------------------------------------

noise_params <- tibble(
  antigen_iso = "IgG",
  nu = 0.5,          # Biological noise parameter
  eps = 0.25,        # Measurement error rate
  y.low = 25,        # Lower detection limit
  y.high = 200000    # Upper detection limit
)

cat("\nNoise parameters:\n")
print(noise_params)

# -----------------------------------------------------------------------------
# 3. Run Small-Scale Simulation (for testing)
# -----------------------------------------------------------------------------
# Run a small number of simulations for demonstration
# In practice, you would run many more replicates (e.g., 1000+)

cat("\nRunning small-scale simulation...\n")

simulation_results <- simulate_seroincidence(
  dmcmc = mock_curve_params,
  nrep = 50,                # Sample size per simulation
  n_sim = 10,               # Number of simulation replicates
  observed = 0.15,          # True incidence rate (per person-year)
  range = c(0.5, 60),       # Age range (years)
  batch_size = 5,           # Simulations per batch
  parallel = FALSE,         # Use FALSE for testing, TRUE for production
  antibodies = c("IgG"),
  cond = noise_params
)

cat("Simulation complete. Number of replicates:", length(simulation_results), "\n")

# -----------------------------------------------------------------------------
# 4. Extract and Summarize Results
# -----------------------------------------------------------------------------

# Extract first simulation result
first_result <- simulation_results[[1]]

cat("\nFirst simulation result structure:\n")
cat("Components:", names(first_result), "\n")

# Extract incidence estimates from all simulations
incidence_estimates <- sapply(simulation_results, function(x) {
  summary(x$est1)$incidence.rate
})

cat("\nIncidence rate estimates across simulations:\n")
cat("True value:", 0.15, "\n")
cat("Mean estimate:", mean(incidence_estimates), "\n")
cat("SD of estimates:", sd(incidence_estimates), "\n")
cat("Range:", range(incidence_estimates), "\n")

# -----------------------------------------------------------------------------
# 5. Generate Summary Table
# -----------------------------------------------------------------------------

# Create final table with confidence intervals
summary_table <- generate_final_table(
  results_list = simulation_results,
  sample_size = 50
)

cat("\nSummary table (first 10 rows):\n")
print(head(summary_table, 10))

# Calculate coverage probability
coverage <- mean(
  summary_table$CI.lwr <= 0.15 & summary_table$CI.upr >= 0.15
)

cat("\nCoverage probability (95% CI contains true value):", 
    round(coverage * 100, 1), "%\n")

# -----------------------------------------------------------------------------
# 6. Create Noise Parameters for Different Regions
# -----------------------------------------------------------------------------

noise_usa <- create_noise_df("MA USA")
noise_ghana <- create_noise_df("Ghana")
noise_niger <- create_noise_df("Niger")

cat("\nNoise parameters for USA:\n")
print(noise_usa)

# -----------------------------------------------------------------------------
# Notes for Production Use:
# -----------------------------------------------------------------------------
# 1. Use real curve parameters from MCMC fitting of longitudinal data
# 2. Increase n_sim to 1000+ for stable estimates
# 3. Use parallel = TRUE for faster computation
# 4. Vary sample sizes (nrep) to assess power
# 5. Compare results across different geographic regions
