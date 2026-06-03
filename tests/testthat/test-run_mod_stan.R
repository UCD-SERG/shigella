test_that("run_mod_stan has required arguments", {
  fn_args <- names(formals(run_mod_stan))
  expect_true("data" %in% fn_args)
  expect_true("model" %in% fn_args)
  expect_true("seed" %in% fn_args)
})

test_that("run_mod_stan completes a minimal fit (slow)", {
  skip_if(
    Sys.getenv("RUN_STAN_TESTS") != "true",
    "Stan tests are skipped unless RUN_STAN_TESTS=true."
  )
  skip_if_not_installed("cmdstanr")

  sim <- sim_correlated_case_data(n = 3, seed = 2026)

  warnings_seen <- character()
  fit <- withCallingHandlers(
    run_mod_stan(
      data          = sim,
      model         = "model_2",
      chains        = 1,
      iter_warmup   = 200,
      iter_sampling = 100,
      adapt_delta   = 0.99,
      max_treedepth = 15,
      refresh       = 0,
      show_messages = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Low-iteration smoke fit may emit convergence warnings — that is expected.
  # Assert that the function ran and returned a valid object; do not assert
  # on the absence of warnings since they are diagnostic signals, not errors.
  expect_s3_class(fit, "sr_model")
})
