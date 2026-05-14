## Example: prep_data_stan()
##
## Convert long-format case data into the structured list that the
## Stan model expects.

set.seed(2026)

sim_data <- sim_correlated_case_data(
  n                 = 5,
  Omega_B           = diag(2),
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

stan_data <- prep_data_stan(sim_data)

str(stan_data)
