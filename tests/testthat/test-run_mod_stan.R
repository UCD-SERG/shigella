test_that("run_mod_stan is callable with expected arguments", {
  expect_true(is.function(run_mod_stan))
  
  fn_args <- names(formals(run_mod_stan))
  expect_true(any(c("data", "case_data") %in% fn_args))
  expect_true(any(c("model", "file_mod") %in% fn_args))
})

test_that("run_mod_stan completes a minimal fit (slow)", {
  skip_if(
    Sys.getenv("RUN_STAN_TESTS") != "true",
    "Stan tests are skipped unless RUN_STAN_TESTS=true."
  )
  skip_if_not_installed("cmdstanr")

  sim <- sim_correlated_case_data(n = 3, seed = 2026)

  fit <- suppressWarnings(run_mod_stan(
    data          = sim,
    model         = "model_2",
    chains        = 1,
    iter_warmup   = 200,
    iter_sampling = 100,
    adapt_delta   = 0.99,
    max_treedepth = 15,
    refresh       = 0,
    show_messages = FALSE
  ))

  expect_s3_class(fit, "sr_model")
})

test_that("run_mod_stan has expected function signature", {
  expect_equal(
    names(formals(run_mod_stan)),
    c("data", "model", "chains", "iter_sampling", "iter_warmup",
      "adapt_delta", "max_treedepth", "seed", "strat", "parallel_chains",
      "with_post", "stan_dir", "compile_dir", "init", "refresh",
      "show_messages", "...")
  )
})
