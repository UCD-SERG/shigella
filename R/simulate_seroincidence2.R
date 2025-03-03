#' Description of the function here.
#' @param nrep Number of repetitions.
#' @param n_sim Number of simulations.
#' @param observed Observed incidence rate.
#' @param range Range for simulation.
#' @return A list of simulated seroincidence results.
#' @export
# Define the simulation function
# Define the simulation function
simulate_seroincidence2 <- function(nrep, n_sim, observed, range = NULL) {
  # Set parallel plan inside function to avoid issues with distributed nodes
  plan(multicore)  # Use multiple cores for parallel processing (works best on HPC)
  
  # Parameters
  dmcmc <- curve_params_shigella_ipab  # Curve parameters
  antibodies <- c("IgG")  # Antigen-isotypes
  lambda <- observed  # Simulated incidence rate per person-year
  
  # Biologic noise distribution
  dlims <- rbind("IgG" = c(min = 0, max = 0.5))
  
  # Noise parameters
  cond <- tibble(
    antigen_iso = c("IgG"),
    nu = c(0.5),  # Biologic noise (nu)
    eps = c(0.25),  # Measurement noise (eps)
    y.low = c(25),  # Low cutoff (llod)
    y.high = c(200000)  # High cutoff (y.high)
  )
  
  # Perform simulations in parallel
  results <- future_map(1:n_sim, function(i) {
    tryCatch({
      # Generate cross-sectional data
      csdata <- sim_pop_data(
        curve_params = dmcmc,
        lambda = lambda,
        n.smpl = nrep,
        age_range = range,
        antigen_isos = antibodies,
        n.mc = 0,
        renew_params = TRUE,  # Use different parameters for each simulation
        add.noise = TRUE,
        noise_limits = dlims,
        format = "long"
      )
      
      # Estimate seroincidence
      est <- est.incidence(
        pop_data = csdata,
        curve_params = dmcmc,
        noise_params = cond,
        lambda_start = 0.1,
        build_graph = TRUE,
        verbose = FALSE,
        print_graph = FALSE,
        antigen_isos = antibodies
      )
      
      # Return results for this simulation
      list(csdata = csdata, est1 = est)
    }, error = function(e) {
      return(list(error = e$message))  # Capture and store errors instead of stopping execution
    })
  }, .options = furrr_options(seed = TRUE))
  
  # Ensure sequential processing after function execution
  plan(sequential)
  results <- results |> 
     structure(
         sample_size = nrep,
         age_range = paste(range, collapse = " - ")
      )
  return(results)
}
