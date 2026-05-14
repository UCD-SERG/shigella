test_that("postprocess_stan_output is callable", {
  expect_true(is.function(postprocess_stan_output))
})

test_that("postprocess_stan_output produces sr_model output (slow)", {
  skip_on_ci()
  skip_if_not_installed("cmdstanr")

  sim <- sim_correlated_case_data(n = 3, seed = 2026)

  ## run_mod_stan() internally calls postprocess_stan_output(),
  ## so verifying its output is verifying the post-processing step.
  fit <- run_mod_stan(
    data          = sim,
    model         = "model_2",
    chains        = 1,
    iter_warmup   = 100,
    iter_sampling = 100,
    refresh       = 0,
    show_messages = FALSE
  )

  expect_s3_class(fit, "sr_model")
  expect_true("priors" %in% names(attributes(fit)))
  expect_true("fitted_residuals" %in% names(attributes(fit)))
})
