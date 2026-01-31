#' Simulate Seroincidence with Biological and Measurement Noise
#'
#' Performs Monte Carlo simulations to estimate seroincidence rates using
#' serodynamics models. Supports parallel processing for efficient computation
#' of multiple simulation replicates.
#'
#' @param dmcmc A data frame or tibble containing antibody decay curve parameters,
#'   typically obtained from MCMC fitting of longitudinal seroresponse data.
#' @param nrep Integer. Number of individuals to simulate in each cross-sectional
#'   sample.
#' @param n_sim Integer. Total number of simulation replicates to perform.
#' @param observed Numeric. The observed incidence rate (per person-year) to use
#'   as the true lambda parameter in simulations.
#' @param range Numeric vector of length 2 specifying the age range for simulated
#'   individuals (e.g., `c(0.5, 60)`). If `NULL`, uses default age distribution.
#' @param batch_size Integer. Number of simulations to run in each parallel batch.
#'   Default is 40.
#' @param parallel Logical. If `TRUE` (default), uses parallel processing via
#'   `future` and `furrr`. If `FALSE`, runs sequentially.
#' @param antibodies Character vector of antigen-isotype combinations to simulate
#'   (e.g., `c("IgG", "IgA")`). Default is `c("IgG")`.
#' @param dlims A named matrix or data frame specifying biological noise limits
#'   for each antigen-isotype. Should have columns `min` and `max`, with row names
#'   matching `antibodies`. Default is `rbind("IgG" = c(min = 0, max = 0.5))`.
#' @param cond A tibble containing noise parameters for each antigen-isotype with
#'   columns: `antigen_iso`, `nu` (decay rate), `eps` (error rate), `y.low`
#'   (lower measurement bound), and `y.high` (upper measurement bound).
#'
#' @return A list of length `n_sim`, where each element contains:
#'   \item{csdata}{Simulated cross-sectional antibody data}
#'   \item{est1}{Estimated seroincidence from the simulated data}
#'
#' @examples
#' \dontrun{
#' # Create mock curve parameters (placeholder)
#' mock_dmcmc <- tibble::tibble(
#'   antigen_iso = rep("IgG", 100),
#'   alpha = rnorm(100, 3, 0.5),
#'   beta = rnorm(100, 0.1, 0.02),
#'   r = rnorm(100, 0.01, 0.002)
#' )
#'
#' # Define noise parameters
#' noise_params <- tibble::tibble(
#'   antigen_iso = "IgG",
#'   nu = 0.5,
#'   eps = 0.25,
#'   y.low = 25,
#'   y.high = 200000
#' )
#'
#' # Run simulations
#' results <- simulate_seroincidence(
#'   dmcmc = mock_dmcmc,
#'   nrep = 200,
#'   n_sim = 100,
#'   observed = 0.1,
#'   range = c(0.5, 60),
#'   batch_size = 20,
#'   parallel = TRUE,
#'   antibodies = c("IgG"),
#'   cond = noise_params
#' )
#' }
#'
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom serodynamics sim_pop_data est.incidence
#' @export
# Define the simulation function
simulate_seroincidence <- function(
    dmcmc, # Curve parameters
    nrep, 
    n_sim, 
    observed, 
    range = NULL, 
    batch_size = 40, parallel = TRUE,
    antibodies = c("IgG"), # Antigen-isotypes
    dlims = rbind("IgG" = c(min = 0, max = 0.5)), # Biologic noise distribution
    # Noise parameters
    cond = tibble(
      antigen_iso = c("IgG"),
      nu = c(0.5),
      eps = c(0.25),
      y.low = c(25),
      y.high = c(200000)
    )
) {
  
  lambda <- observed # Simulated incidence rate per person-year
  
  
  
  # Calculate number of batches
  n_batches <- ceiling(n_sim / batch_size)
  
  # Parallel processing: Use future_map for efficiency
  if (parallel) {
    plan(multisession, workers = max(1, future::availableCores() - 1))
    
    results <- future_map(1:n_batches, function(batch) {
      # Run simulations within each batch
      replicate(batch_size, expr = {
        csdata <- sim_pop_data(
          curve_params = dmcmc,
          lambda = lambda,
          n_samples = nrep,
          age_range = range,
          antigen_isos = antibodies,
          n_mcmc_samples = 0,
          renew_params = FALSE, # Keep the same params for speed
          add_noise = TRUE,
          noise_limits = dlims,
          format = "long"
        )
        
        # Estimate seroincidence
        est <- est.incidence(
          pop_data = csdata,
          curve_params = dmcmc,
          noise_params = cond,
          lambda_start = 0.1,
          build_graph = FALSE, # Disable visualization for speed
          verbose = FALSE,
          print_graph = FALSE,
          antigen_isos = antibodies
        )
        
        list(csdata = csdata, est1 = est)
      }, simplify = FALSE) # Keep list structure
    }, .options = furrr_options(seed = TRUE)) # Set seed for reproducibility
    
    # Flatten nested list
    results <- unlist(results, recursive = FALSE)
  } else {
    # Sequential execution (for debugging or single-core use)
    results <- lapply(1:n_sim, function(i) {
      csdata <- sim.cs(
        curve_params = dmcmc,
        lambda = lambda,
        n.smpl = nrep,
        age_range = range,
        antigen_isos = antibodies,
        n.mc = 0,
        renew_params = FALSE,
        add.noise = TRUE,
        noise_limits = dlims,
        format = "long"
      )
      
      est <- est.incidence(
        pop_data = csdata,
        curve_params = dmcmc,
        noise_params = cond,
        lambda_start = 0.1,
        build_graph = FALSE,
        verbose = FALSE,
        print_graph = FALSE,
        antigen_isos = antibodies
      )
      
      list(csdata = csdata, est1 = est)
    })
  }
  
  return(results)
}