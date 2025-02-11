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