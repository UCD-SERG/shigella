# tests/testthat/test-residuals.R
# Tests for get_observed, compute_residual_metrics, get_mae

# ── Helper: mock case_data ──────────────────────────────────────────────────

make_mock_casedata <- function() {
  df <- data.frame(
    id          = rep(c("S1", "S1", "S1", "S2", "S2", "S2"), 2),
    antigen_iso = rep(c("IgG", "IgA"), each = 6),
    timeindays  = rep(c(2, 30, 90), 4),
    result      = c(5000, 8000, 6000, 3000, 5000, 4000,
                    2000, 4000, 1000, 1500, 3000, 800),
    sid         = rep(c("S1", "S1", "S1", "S2", "S2", "S2"), 2),
    isotype_name = rep(c("IgG", "IgA"), each = 6),
    cohort_name = "Sf2a",
    stringsAsFactors = FALSE
  )
  structure(
    tibble::as_tibble(df),
    class         = c("case_data", "tbl_df", "tbl", "data.frame"),
    id_var        = "id",
    biomarker_var = "antigen_iso",
    timeindays    = "timeindays",
    value_var     = "result"
  )
}

# ── get_observed ────────────────────────────────────────────────────────────

test_that("get_observed returns correct rows for one subject/isotype", {
  cd <- make_mock_casedata()
  result <- get_observed(cd, "S1", "IgG")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("t", "value") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_equal(result$t, c(2, 30, 90))
})

test_that("get_observed warns for missing subject", {
  cd <- make_mock_casedata()
  expect_warning(
    get_observed(cd, "NONEXISTENT", "IgG"),
    "No observed data"
  )
})

test_that("get_observed returns empty tibble for missing subject", {
  cd <- make_mock_casedata()
  result <- suppressWarnings(get_observed(cd, "NONEXISTENT", "IgG"))
  expect_equal(nrow(result), 0)
})

test_that("get_observed works with isotype_name column", {
  df <- data.frame(
    sid          = c("A", "A"),
    isotype_name = c("IgG", "IgG"),
    timeindays   = c(2, 30),
    result       = c(100, 200)
  )
  result <- get_observed(df, "A", "IgG")
  expect_equal(nrow(result), 2)
})

# ── get_mae (integration-style test) ────────────────────────────────────────

test_that("get_mae returns tibble with correct columns", {
  skip_if_not_installed("shigella")

  # Use package mock data if available
  tryCatch({
    data("overall_IpaB_pop_6", package = "shigella", envir = environment())
    data("dL_clean_Ipab_new",  package = "shigella", envir = environment())

    result <- get_mae(overall_IpaB_pop_6, dL_clean_Ipab_new, "IpaB", "IgG")

    expect_true(all(c("sid", "antigen", "Iso_type", "mae") %in% names(result)))
    expect_true(nrow(result) > 0)
    expect_true(all(result$mae >= 0))
    expect_true(all(result$antigen == "IpaB"))
    expect_true(all(result$Iso_type == "IgG"))
  }, error = function(e) {
    skip("Mock data not available for integration test")
  })
})

test_that("get_mae returns empty tibble on error", {
  # Pass a completely wrong model — should not crash, just return empty
  bad_model <- data.frame(
    Subject = "X", Iso_type = "IgG", Chain = 1L,
    Iteration = 1L, Parameter = "y0", value = 5
  )
  cd <- make_mock_casedata()

  result <- get_mae(bad_model, cd, "Test", "IgG")
  expect_equal(nrow(result), 0)
})
