## Example: run_mod_stan()
##
## Fit the Chapter 2 Kronecker Stan model on a small synthetic
## dataset. Uses minimal MCMC settings so the example completes
## quickly. For realistic settings, see the Phase 2 simulation
## scripts (run on Shiva HPC).

if (requireNamespace("cmdstanr", quietly = TRUE)) {

  set.seed(2026)

  sim_data <- sim_correlated_case_data(
    n                 = 5,
    Omega_B           = matrix(c(1, 0.6, 0.6, 1), nrow = 2),
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )

  fit <- run_mod_stan(
    data          = sim_data,
    model         = "model_2",
    chains        = 1,
    iter_warmup   = 200,
    iter_sampling = 200,
    refresh       = 0,
    show_messages = FALSE
  )

  class(fit)
  head(fit)
}
