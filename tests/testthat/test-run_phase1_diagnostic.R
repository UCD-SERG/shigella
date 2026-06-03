test_that("run_phase1_diagnostic has required arguments", {
  fn_args <- names(formals(shigella:::run_phase1_diagnostic))
  expect_true("n" %in% fn_args)
  expect_true("tag" %in% fn_args)
  expect_true("compile_dir" %in% fn_args)
  expect_true("phase0_dir" %in% fn_args)
})

test_that("run_phase1_diagnostic runs with minimal Stan settings", {
  skip_if_not(Sys.getenv("RUN_STAN_TESTS") == "true",
              "Skipping Stan test (set RUN_STAN_TESTS=true to enable)")

  out_dir  <- file.path(tempdir(), "phase1_test")
  ph0_dir  <- file.path(tempdir(), "phase0_test")
  on.exit(unlink(c(out_dir, ph0_dir), recursive = TRUE))

  result <- run_phase1_diagnostic(
    n             = 3,
    iter_warmup   = 50,
    iter_sampling = 50,
    tag           = "test",
    output_dir    = out_dir,
    phase0_dir    = ph0_dir,
    chains        = 1L
  )

  expect_true(is.list(result) || is.null(result))
  if (!is.null(result)) {
    expect_equal(result$status, "OK")
    expect_equal(result$n_subjects, 3L)
    # run_phase1_diagnostic saves as one_fit_{tag}_jobid_{timestamp}.rds
    rds_files <- list.files(out_dir, pattern = "^one_fit_test.*\\.rds$")
    expect_true(length(rds_files) > 0)
  }
})
