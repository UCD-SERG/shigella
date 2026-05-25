test_that("postprocess_stan_output is callable", {
  expect_true(is.function(postprocess_stan_output))
})

test_that("summarize_matrix_array returns a list of matrices", {
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
  result <- shigella:::.summarize_matrix_array(arr_data, "Omega_P",
                                              n_arr = 2L, n_row = 3L, n_col = 3L)

  expect_type(result, "list")
  expect_length(result, 2L)
  expect_equal(dim(result[[1]]), c(3L, 3L))
  expect_equal(dim(result[[2]]), c(3L, 3L))
  # Each cell should equal the draw value (only one draw, so median = the value)
  expect_equal(result[[1]][1, 1],
               arr_data[1, 1, "Omega_P[1,1,1]"],
               tolerance = 1e-10)
})

# ─── Mock Stan fit for fast (non-Stan) pipeline tests ───────────────────────

# Builds a minimal mock CmdStanMCMC-like list with just enough interface to
# satisfy postprocess_stan_output: $draws(variables) and $metadata().
# Avoids running Stan in routine devtools::test() calls.
.make_mock_stan_fit <- function(N = 2L, K = 2L,
                                param_names = c("y0", "y1", "t1",
                                                "alpha", "shape"),
                                n_iter = 4L, n_chains = 2L) {
  P <- length(param_names)

  # Parameter variables: pname[subj, antigen] for all combinations
  param_vars <- unlist(lapply(param_names, function(p) {
    sprintf("%s[%d,%d]", p,
            rep(seq_len(N), each = K),
            rep(seq_len(K), N))
  }))

  # model_1 Omega_P: array[K] corr_matrix[P] -> Omega_P[k,i,j]
  omega_P_vars <- unlist(lapply(seq_len(K), function(k) {
    as.vector(outer(seq_len(P), seq_len(P),
                    function(i, j) sprintf("Omega_P[%d,%d,%d]", k, i, j)))
  }))

  log_lik_vars <- sprintf("log_lik[%d]", seq_len(N))

  all_vars <- c(param_vars, omega_P_vars, log_lik_vars)

  set.seed(42)
  raw_arr <- array(
    abs(rnorm(n_iter * n_chains * length(all_vars), mean = 1, sd = 0.2)),
    dim = c(n_iter, n_chains, length(all_vars))
  )
  dimnames(raw_arr) <- list(NULL, NULL, all_vars)
  draws_full <- posterior::as_draws_array(raw_arr)

  list(
    draws = function(variables = NULL, ...) {
      if (is.null(variables)) return(draws_full)
      all_var_names <- posterior::variables(draws_full)
      matched <- unlist(lapply(variables, function(v) {
        grep(paste0("^", v, "(\\[|$)"), all_var_names, value = TRUE)
      }))
      if (length(matched) == 0) {
        stop("No variables matched: ", paste(variables, collapse = ", "))
      }
      posterior::subset_draws(draws_full, variable = matched)
    },
    metadata = function() {
      list(stan_variables = c(param_names, "Omega_P", "log_lik"))
    }
  )
}

test_that("postprocess_stan_output processes mock draws without Stan (model_1)", {
  skip_if_not_installed("posterior")

  ids      <- c("s1", "s2")
  antigens <- c("IgG", "IgA")
  param_names <- c("y0", "y1", "t1", "alpha", "shape")
  n_iter   <- 4L
  n_chains <- 2L

  mock_fit <- .make_mock_stan_fit(
    N           = length(ids),
    K           = length(antigens),
    param_names = param_names,
    n_iter      = n_iter,
    n_chains    = n_chains
  )

  result <- postprocess_stan_output(
    stan_fit       = mock_fit,
    ids            = ids,
    antigens       = antigens,
    model          = "model_1",
    stratification = "test_stratum"
  )

  # Returns named list with sr_tibble and cov_summaries
  expect_type(result, "list")
  expect_named(result, c("sr_tibble", "cov_summaries"))

  # sr_tibble has expected columns
  expect_s3_class(result$sr_tibble, "tbl_df")
  expected_cols <- c("Iteration", "Chain", "Parameter", "Iso_type",
                     "Stratification", "Subject", "value")
  expect_true(all(expected_cols %in% names(result$sr_tibble)))

  # All 5 parameters present for both subjects and antigens
  expect_equal(length(unique(result$sr_tibble$Parameter)), length(param_names))
  expect_equal(sort(unique(result$sr_tibble$Subject)), sort(ids))
  expect_equal(sort(unique(result$sr_tibble$Iso_type)), sort(antigens))

  # Row count: n_params × N × K × (n_iter × n_chains)
  expect_equal(
    nrow(result$sr_tibble),
    length(param_names) * length(ids) * length(antigens) * n_iter * n_chains
  )

  # Numerical sanity: all values finite, ESS-free check
  expect_true(all(is.finite(result$sr_tibble$value)))

  # Omega_P extracted (model_1 produces a named list of K matrices)
  expect_true("Omega_P" %in% names(result$cov_summaries))
  omega_P <- result$cov_summaries$Omega_P
  expect_type(omega_P, "list")
  expect_equal(length(omega_P), length(antigens))
  for (mat in omega_P) {
    expect_equal(dim(mat), c(length(param_names), length(param_names)))
    expect_true(all(is.finite(mat)))
  }
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
