## Example: prep_data_stan()
##
## Convert a case_data object directly into the list format that
## Stan models expect. Internally calls serodynamics::prep_data()
## with add_newperson = FALSE.

library(shigella)

sim <- sim_correlated_case_data(n = 5, seed = 2026)

stan_data <- prep_data_stan(sim)

str(stan_data)
