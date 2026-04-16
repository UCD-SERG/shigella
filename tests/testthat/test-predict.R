# tests/testthat/test-predict.R
# Tests for predict_posterior_at_times, get_prediction_summary

# ── Helper: minimal model with individual draws ─────────────────────────────

make_indiv_model <- function(sids = c("S1", "S2"),
                             n_iter = 5L,
                             n_chains = 1L) {
  params <- c("y0", "y1", "t1", "alpha", "shape")
  # Realistic log-scale values for IpaB-like kinetics
  base_vals <- c(y0 = 5.0, y1 = 10.0, t1 = -0.5, alpha = -15.0, shape = 0.3)

  rows <- list()
  for (sid in sids) {
    for (iso in c("IgG", "IgA")) {
      for (ch in seq_len(n_chains)) {
        for (it in seq_len(n_iter)) {
          vals <- base_vals + stats::rnorm(5, 0, 0.1)
          for (p in seq_along(params)) {
            rows[[length(rows) + 1L]] <- data.frame(
              Iteration      = it,
              Chain          = ch,
              Parameter      = params[p],
              Iso_type       = iso,
              Stratification = "None",
              Subject        = sid,
              value          = vals[p],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  structure(
    dplyr::bind_rows(rows),
    class = c("sr_model", "tbl_df", "tbl", "data.frame")
  )
}

# ── predict_posterior_at_times ───────────────────────────────────────────────

test_that("predict_posterior_at_times returns a tibble with res column", {
  model <- make_indiv_model("S1", n_iter = 3L)
  result <- predict_posterior_at_times(
    model = model, ids = "S1", antigen_iso = "IgG",
    times = c(0, 30, 90)
  )
  expect_s3_class(result, "data.frame")
  expect_true("res" %in% names(result))
  expect_true("t" %in% names(result))
  expect_true(all(result$res > 0))
})

test_that("predict_posterior_at_times has correct dimensions", {
  model <- make_indiv_model("S1", n_iter = 4L, n_chains = 1L)
  times <- c(0, 50, 100)
  result <- predict_posterior_at_times(
    model = model, ids = "S1", antigen_iso = "IgG", times = times
  )
  # 4 iterations × 1 chain × 3 timepoints = 12 rows
  expect_equal(nrow(result), 4 * 1 * 3)
})

test_that("predict_posterior_at_times works with multiple subjects", {
  model <- make_indiv_model(c("S1", "S2"), n_iter = 3L)
  result <- predict_posterior_at_times(
    model = model, ids = c("S1", "S2"), antigen_iso = "IgG",
    times = c(0, 100)
  )
  subjects_in_result <- unique(as.character(result$id))
  expect_true("S1" %in% subjects_in_result)
  expect_true("S2" %in% subjects_in_result)
})

# ── get_prediction_summary ──────────────────────────────────────────────────

test_that("get_prediction_summary returns med/lo/hi columns", {
  model <- make_indiv_model("S1", n_iter = 10L)
  result <- get_prediction_summary(model, "S1", "IgG",
                                   times = seq(0, 100, by = 50))
  expect_true(all(c("t", "med", "lo", "hi") %in% names(result)))
  expect_equal(nrow(result), 3)  # 3 time points
  # median should be between lo and hi
  expect_true(all(result$med >= result$lo))
  expect_true(all(result$med <= result$hi))
})

test_that("get_prediction_summary warns for missing subject", {
  model <- make_indiv_model("S1")
  expect_warning(
    get_prediction_summary(model, "NONEXISTENT", "IgG", times = c(0)),
    "No data"
  )
})

test_that("get_prediction_summary returns empty tibble for missing subject", {
  model <- make_indiv_model("S1")
  result <- suppressWarnings(
    get_prediction_summary(model, "NONEXISTENT", "IgG", times = c(0))
  )
  expect_equal(nrow(result), 0)
})
