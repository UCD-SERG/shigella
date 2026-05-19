## Example: postprocess_stan_output()
##
## Convert a raw cmdstanr fit object into a list with sr_tibble
## (tidy parameter draws) and cov_summaries (covariance matrices).
## Requires a compiled cmdstan installation.

if (interactive()) {
if (requireNamespace("cmdstanr", quietly = TRUE)) {

  set.seed(2026)

  sim_data <- sim_correlated_case_data(
    n                 = 5,
    omega_B           = diag(2),
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
  processed <- postprocess_stan_output(
    stan_fit  = raw_fit,
    ids       = attr(stan_data, "ids"),
    antigens  = attr(stan_data, "antigens"),
    model     = "model_2"
  )

  class(processed$sr_tibble)
  names(processed$cov_summaries)
}
}
