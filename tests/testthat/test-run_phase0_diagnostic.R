test_that("run_phase0_diagnostic has required arguments", {
  fn_args <- names(formals(shigella:::run_phase0_diagnostic))
  expect_true("n" %in% fn_args)
  expect_true("tag" %in% fn_args)
  expect_true("compile_dir" %in% fn_args)
})

test_that("run_phase0_diagnostic runs with minimal Stan settings", {
  skip_if_not(Sys.getenv("RUN_STAN_TESTS") == "true",
              "Skipping Stan test (set RUN_STAN_TESTS=true to enable)")

  out_dir <- file.path(tempdir(), "phase0_test")
  on.exit(unlink(out_dir, recursive = TRUE))

  result <- run_phase0_diagnostic(
    n             = 3,
    iter_warmup   = 50,
    iter_sampling = 50,
    tag           = "test",
    output_dir    = out_dir,
    chains        = 1L
  )

  expect_true(is.list(result) || is.null(result))
  if (!is.null(result)) {
    expect_equal(result$status, "OK")
    expect_equal(result$n_subjects, 3L)
    expect_true(file.exists(file.path(out_dir, "one_fit_test.rds")))
  }
})
