test_that("postprocess_stan_output is callable", {
  expect_true(is.function(postprocess_stan_output))
})

test_that("summarize_array_of_matrix_draws returns a list of matrices", {
  # Build a minimal mock draws_array for array[2] corr_matrix[3] Omega_P
  # Variable names: Omega_P[k,i,j] for k in 1:2, i in 1:3, j in 1:3
  var_names <- c()
  for (k in 1:2) {
    for (i in 1:3) {
      for (j in 1:3) {
        var_names <- c(var_names, sprintf("Omega_P[%d,%d,%d]", k, i, j))
      }
    }
  }
  set.seed(1)
  arr_data <- array(
    runif(1 * 1 * length(var_names)),
    dim = c(1L, 1L, length(var_names)),
    dimnames = list(
      iteration = "1",
      chain     = "1",
      variable  = var_names
    )
  )
  result <- summarize_array_of_matrix_draws(arr_data, "Omega_P",
                                            n_arr = 2L, nrow = 3L, ncol = 3L)

  expect_type(result, "list")
  expect_length(result, 2L)
  expect_equal(dim(result[[1]]), c(3L, 3L))
  expect_equal(dim(result[[2]]), c(3L, 3L))
  # Each cell should equal the draw value (only one draw, so median = the value)
  expect_equal(result[[1]][1, 1],
               arr_data[1, 1, "Omega_P[1,1,1]"],
               tolerance = 1e-10)
})

test_that("postprocess_stan_output produces sr_model output — model_2 (slow)", {
  skip_if(
    Sys.getenv("RUN_STAN_TESTS") != "true",
    "Stan tests are skipped unless RUN_STAN_TESTS=true."
  )
  skip_if_not_installed("cmdstanr")

  sim <- sim_correlated_case_data(n = 3, seed = 2026)

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

test_that("postprocess_stan_output extracts per-biomarker Omega_P — model_1 (slow)", {
  skip_if(
    Sys.getenv("RUN_STAN_TESTS") != "true",
    "Stan tests are skipped unless RUN_STAN_TESTS=true."
  )
  skip_if_not_installed("cmdstanr")

  sim <- sim_correlated_case_data(n = 3, seed = 2026)

  fit <- run_mod_stan(
    data          = sim,
    model         = "model_1",
    chains        = 1,
    iter_warmup   = 100,
    iter_sampling = 100,
    refresh       = 0,
    show_messages = FALSE
  )

  expect_s3_class(fit, "sr_model")

  # run_mod_stan flattens cov_summaries into individual top-level attributes
  expect_true("Omega_P" %in% names(attributes(fit)))
  omega_P <- attr(fit, "Omega_P")

  # model_1 Omega_P should be a named list (one 5x5 matrix per biomarker)
  expect_type(omega_P, "list")
  expect_equal(length(omega_P), 2L)  # K=2 biomarkers (biomarker_1, biomarker_2)
  expect_equal(names(omega_P), c("biomarker_1", "biomarker_2"))
  for (mat in omega_P) {
    expect_equal(dim(mat), c(5L, 5L))
    expect_equal(rownames(mat), c("y0", "y1", "t1", "alpha", "shape"))
    expect_equal(colnames(mat), c("y0", "y1", "t1", "alpha", "shape"))
  }

  # model_1 should NOT have Kronecker-only summaries
  expect_false("Omega_B" %in% names(attributes(fit)))
  expect_false("Omega_eps" %in% names(attributes(fit)))
})
