## Code to prepare mock datasets for testing and examples
## This script generates mock data that mimics the structure expected by shigella functions

set.seed(2024)

# Mock posterior draws dataset -------------------------------------------------
# Structure: Subject, Iso_type, Chain, Iteration, Parameter, value
# Parameters: y0, y1, t1, alpha, shape (rho)

n_chains <- 3
n_iter <- 100
subjects <- c("newperson", "SOSAR-22008", "SOSAR-22015")
isotypes <- c("IgG", "IgA")
parameters <- c("y0", "y1", "t1", "alpha", "shape")

mock_posterior_draws <- expand.grid(
  Subject = subjects,
  Iso_type = isotypes,
  Chain = 1:n_chains,
  Iteration = 1:n_iter,
  Parameter = parameters,
  stringsAsFactors = FALSE
)

# Generate plausible parameter values
mock_posterior_draws$value <- NA
for (i in seq_len(nrow(mock_posterior_draws))) {
  param <- mock_posterior_draws$Parameter[i]
  mock_posterior_draws$value[i] <- switch(param,
    y0 = runif(1, 100, 500),      # baseline antibody level
    y1 = runif(1, 1000, 5000),    # peak antibody level
    t1 = runif(1, 5, 20),         # time to peak (days)
    alpha = runif(1, 0.001, 0.05), # decay rate
    shape = runif(1, 0.5, 2.0)    # decay shape (rho)
  )
}

# Mock case dataset -------------------------------------------------------------
# Structure compatible with serodynamics::as_case_data()
# Required: id, antigen_iso, time (in days), value columns

n_timepoints <- 6
timepoints <- c(0, 7, 14, 30, 60, 90)
ids_case <- c("SOSAR-22008", "SOSAR-22015", "SOSAR-22020", "SOSAR-22025")

mock_case_data_raw <- expand.grid(
  id = ids_case,
  antigen_iso = c("IgG", "IgA"),
  timepoint = timepoints,
  stringsAsFactors = FALSE
)

# Generate mock antibody values using a simple kinetics model
mock_case_data_raw$value <- apply(mock_case_data_raw, 1, function(row) {
  t <- as.numeric(row["timepoint"])
  # Simple rise-decay pattern with noise
  y0 <- runif(1, 100, 300)
  y1 <- runif(1, 1000, 3000)
  t1 <- 10
  alpha <- 0.02
  shape <- 1.2
  
  if (t <= t1) {
    pred <- y0 + (y1 - y0) * (t / t1)
  } else {
    pred <- y0 + (y1 - y0) * exp(-alpha * ((t - t1)^shape))
  }
  
  # Add measurement noise (20% CV)
  pred * rnorm(1, 1, 0.2)
})

# Convert to serodynamics case_data format
mock_case_data <- mock_case_data_raw
attr(mock_case_data, "timeindays") <- "timepoint"
attr(mock_case_data, "value_var") <- "value"
class(mock_case_data) <- c("case_data", "data.frame")

# Save datasets -----------------------------------------------------------------
if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}
save(mock_posterior_draws, file = file.path("data", "mock_posterior_draws.rda"))
save(mock_case_data, file = file.path("data", "mock_case_data.rda"))

message("Mock datasets created successfully!")
message("- mock_posterior_draws: ", nrow(mock_posterior_draws), " rows")
message("- mock_case_data: ", nrow(mock_case_data), " rows")
