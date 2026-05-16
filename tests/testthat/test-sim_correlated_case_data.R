test_that("sim_correlated_case_data returns a case_data object", {
  sim <- sim_correlated_case_data(n = 5, seed = 2026)

  expect_s3_class(sim, "case_data")
  expect_true(nrow(sim) > 0)
  expect_true(all(c("id", "timeindays", "antigen_iso", "value") 
                  %in% names(sim)))
})

test_that("sim_correlated_case_data attaches the truth attributes", {
  omega_B <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sim <- sim_correlated_case_data(n = 5, omega_B = omega_B, seed = 2026)

  truth <- attr(sim, "truth")
  expect_type(truth, "list")
  expect_equal(truth$omega_B, omega_B)
  expect_true(!is.null(attr(sim, "theta_true")))
})

test_that("sim_correlated_case_data supports n = 1", {
  sim <- sim_correlated_case_data(n = 1, seed = 2026)

  expect_s3_class(sim, "case_data")
  expect_equal(length(unique(sim$id)), 1)
})

test_that("sim_correlated_case_data theta_true structure is stable", {
  sim <- sim_correlated_case_data(n = 5, seed = 2026)
  theta <- attr(sim, "theta_true")

  expect_equal(dim(theta), c(5L, 5L, 2L))
  expect_equal(
    dimnames(theta),
    list(
      subject  = as.character(1:5),
      param    = c("log_y0", "log_y1m0", "log_t1", "log_alpha", "log_rm1"),
      biomarker = c("biomarker_1", "biomarker_2")
    )
  )

  truth <- attr(sim, "truth")
  expect_equal(
    names(truth),
    c("mu", "tau_P", "tau_B", "tau_eps", "omega_P", "omega_B", "omega_eps",
      "sigma_P", "sigma_B", "sigma_eps")
  )
})
