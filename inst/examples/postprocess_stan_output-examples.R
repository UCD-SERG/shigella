## Example: postprocess_stan_output()
##
## Convert a raw cmdstanr fit object into the tidy `sr_model` format
## with priors and fitted_residuals attached as attributes.

if (requireNamespace("cmdstanr", quietly = TRUE)) {

  set.seed(2026)

  sim_data <- sim_correlated_case_data(
    n                 = 5,
    Omega_B           = diag(2),
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )

  stan_data <- prep_data_stan(sim_data)
  priors    <- prep_priors_stan(model = "model_2")

  ## Compile and sample directly (without run_mod_stan wrapper)
  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "model_2.stan", package = "shigella")
  )
  raw_fit <- mod$sample(
    data            = c(stan_data, priors),
    chains          = 1,
    iter_warmup     = 200,
    iter_sampling   = 200,
    refresh         = 0,
    show_messages   = FALSE
  )

  ## Post-process the raw fit into tidy sr_model format
  tidy_fit <- postprocess_stan_output(
    stan_fit      = raw_fit,
    original_data = sim_data,
    priors        = priors
  )

  class(tidy_fit)
  head(tidy_fit)
}
