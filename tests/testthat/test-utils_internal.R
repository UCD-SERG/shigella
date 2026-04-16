# tests/testthat/test-utils_internal.R
# Tests for .ab(), .bt(), .as_case_data(), and helper functions

# ── .ab() — Two-phase antibody trajectory (LOG-scale inputs) ────────────────

test_that(".ab() returns baseline at t = 0", {
 # log-scale inputs: y0=log(200), y1=log(1800), t1=log(10), alpha=log(0.02), shape=log(0.2)
 result <- shigella:::.ab(
    t = 0,
    y0 = log(200),
    y1 = log(1800),    # peak = exp(y0) + exp(y1) = 200 + 1800 = 2000
    t1 = log(10),
    alpha = log(0.02),
    shape = log(0.2)   # shape_nat = exp(log(0.2)) + 1 = 1.2
  )
  expect_equal(result, 200, tolerance = 1e-6)
})

test_that(".ab() returns peak at t = t1", {
  result <- shigella:::.ab(
    t = 10,
    y0 = log(200),
    y1 = log(1800),
    t1 = log(10),
    alpha = log(0.02),
    shape = log(0.2)
  )
  # peak = exp(log(200)) + exp(log(1800)) = 200 + 1800 = 2000
  expect_equal(result, 2000, tolerance = 1e-6)
})

test_that(".ab() shows exponential rise for t < t1", {
  y0_log <- log(200)
  y1_log <- log(1800)
  t1_log <- log(10)

  result <- shigella:::.ab(
    t = 5, y0 = y0_log, y1 = y1_log, t1 = t1_log,
    alpha = log(0.02), shape = log(0.2)
  )

  # Expected: y0_nat * exp(mu_y * t) where mu_y = log(peak/y0_nat) / t1_nat
  y0_nat <- 200
  peak   <- 2000
  t1_nat <- 10
  mu_y   <- log(peak / y0_nat) / t1_nat
  expected <- y0_nat * exp(mu_y * 5)

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that(".ab() decays after peak (t > t1)", {
  result_at_peak <- shigella:::.ab(
    t = 10, y0 = log(200), y1 = log(1800), t1 = log(10),
    alpha = log(0.02), shape = log(0.2)
  )
  result_after <- shigella:::.ab(
    t = 100, y0 = log(200), y1 = log(1800), t1 = log(10),
    alpha = log(0.02), shape = log(0.2)
  )
  expect_lt(result_after, result_at_peak)
})

test_that(".ab() is vectorized over t", {
  times <- c(0, 5, 10, 50, 100)
  result <- shigella:::.ab(
    t = times, y0 = log(200), y1 = log(1800), t1 = log(10),
    alpha = log(0.02), shape = log(0.2)
  )
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))
})

test_that(".ab() is vectorized over parameters", {
  y0_vec <- log(c(200, 300))
  y1_vec <- log(c(1800, 2700))
  result <- shigella:::.ab(
    t = 5, y0 = y0_vec, y1 = y1_vec,
    t1 = log(c(10, 10)), alpha = log(c(0.02, 0.02)),
    shape = log(c(0.2, 0.2))
  )
  expect_length(result, 2)
})

# ── .bt() — Power-law decay kernel ──────────────────────────────────────────

test_that(".bt() returns 1 at t = t1 (no decay yet)", {
  result <- shigella:::.bt(t = 10, t1 = 10, alpha = 0.02, shape = 1.2)
  expect_equal(result, 1.0, tolerance = 1e-10)
})

test_that(".bt() returns value < 1 for t > t1", {
  result <- shigella:::.bt(t = 100, t1 = 10, alpha = 0.02, shape = 1.2)
  expect_lt(result, 1.0)
  expect_gt(result, 0.0)
})

# ── .as_case_data() — Case data constructor ─────────────────────────────────

test_that(".as_case_data() sets correct class and attributes", {
  df <- data.frame(
    my_id  = c("A", "A", "B"),
    my_iso = c("IgG", "IgG", "IgA"),
    my_t   = c(2, 30, 7),
    my_val = c(100, 200, 150)
  )
  result <- shigella:::.as_case_data(
    df,
    id_var = "my_id", biomarker_var = "my_iso",
    time_in_days = "my_t", value_var = "my_val"
  )
  expect_s3_class(result, "case_data")
  expect_equal(attr(result, "id_var"), "my_id")
  expect_equal(attr(result, "value_var"), "my_val")
  expect_true("id" %in% names(result))
})

# ── .get_timeindays_var / .get_values_var ───────────────────────────────────

test_that(".get_timeindays_var reads attribute", {
  df <- structure(
    data.frame(timeindays = 1:3),
    timeindays = "timeindays"
  )
  expect_equal(shigella:::.get_timeindays_var(df), "timeindays")
})

test_that(".get_values_var reads attribute", {
  df <- structure(
    data.frame(result = 1:3),
    value_var = "result"
  )
  expect_equal(shigella:::.get_values_var(df), "result")
})
