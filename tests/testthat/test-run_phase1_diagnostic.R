test_that("run_phase1_diagnostic has the expected function signature", {
  args <- names(formals(run_phase1_diagnostic))
  expect_equal(args, c(
    "n", "iter_warmup", "iter_sampling", "tag",
    "output_dir", "phase0_dir", "true_rho_B", "seed",
    "chains", "adapt_delta", "max_treedepth"
  ))
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
    expect_true(file.exists(file.path(out_dir, "one_fit_test.rds")))
  }
})
