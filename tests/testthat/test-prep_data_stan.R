test_that("prep_data_stan accepts a case_data object directly", {
  sim <- sim_correlated_case_data(n = 5, seed = 2026)
  
  stan_data <- prep_data_stan(sim)
  
  expect_type(stan_data, "list")
  expect_true(all(c("N", "K", "P", "max_obs", "n_obs",
                    "time_obs", "log_y") %in% names(stan_data)))
  expect_equal(stan_data$N, 5)  # newperson dropped via add_newperson=FALSE
})

test_that("prep_data_stan accepts a prepped_jags_data object", {
  sim <- sim_correlated_case_data(n = 5, seed = 2026)
  prepped <- serodynamics::prep_data(sim, add_newperson = FALSE)
  
  stan_data <- prep_data_stan(prepped)
  
  expect_type(stan_data, "list")
  expect_equal(stan_data$N, 5)
})

test_that("prep_data_stan rejects invalid input class", {
  expect_error(prep_data_stan(list(a = 1, b = 2)))
})
