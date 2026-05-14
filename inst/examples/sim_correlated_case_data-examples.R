## Example: sim_correlated_case_data()
##
## Generate synthetic Shigella antibody-kinetics data with a known
## IgG-IgA correlation rho_B.

set.seed(2026)

Omega_B <- matrix(c(1.0, 0.6,
                    0.6, 1.0), nrow = 2)

sim_data <- sim_correlated_case_data(
  n                 = 5,
  Omega_B           = Omega_B,
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

head(sim_data)
