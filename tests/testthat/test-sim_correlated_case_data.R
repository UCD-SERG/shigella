test_that("sim_correlated_case_data returns a case_data object", {
  sim <- sim_correlated_case_data(n = 5, seed = 2026)

  expect_s3_class(sim, "case_data")
  expect_true(nrow(sim) > 0)
  expect_true(all(c("id", "timeindays", "antigen_iso", "value") %in% names(sim)))
})

test_that("sim_correlated_case_data attaches the truth attributes", {
  Omega_B <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sim <- sim_correlated_case_data(n = 5, Omega_B = Omega_B, seed = 2026)

  truth <- attr(sim, "truth")
  expect_type(truth, "list")
  expect_equal(truth$Omega_B, Omega_B)
  expect_true(!is.null(attr(sim, "theta_true")))
})
