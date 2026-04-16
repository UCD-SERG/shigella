# tests/testthat/test-kinetics.R
# Tests for prep_pop_params, transform_pop_params, summarise_pop_params,
# fmt_mci, format_param_table

# ── Helper: minimal mock model object ───────────────────────────────────────

make_tiny_model <- function() {
  pop <- tibble::tibble(
    Iteration            = rep(1:4, each = 5 * 2),
    Chain                = rep(1L, 40),
    Parameter            = rep(c("y0", "y1", "t1", "alpha", "shape"), 8),
    Iso_type             = rep(rep(c("IgG", "IgA"), each = 5), 4),
    Stratification       = "None",
    Population_Parameter = "mu.par",
    value                = c(
      # IgG iter 1-4
      5.0, 10.0, -0.5, -15, 0.3,
      5.1, 10.1, -0.4, -14, 0.2,
      4.9,  9.9, -0.6, -16, 0.4,
      5.2, 10.2, -0.3, -13, 0.1,
      # IgA iter 1-4
      3.5,  9.0, -0.8, -12, 0.2,
      3.6,  9.1, -0.7, -11, 0.1,
      3.4,  8.9, -0.9, -13, 0.3,
      3.7,  9.2, -0.6, -10, 0.15
    )
  )
  structure(
    tibble::tibble(),
    class = c("sr_model", "tbl_df", "tbl", "data.frame"),
    population_params = pop
  )
}

# ── prep_pop_params ─────────────────────────────────────────────────────────

test_that("prep_pop_params returns correct columns", {
  model <- make_tiny_model()
  result <- prep_pop_params(model, "IpaB")
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("Iteration", "Chain", "antigen", "Iso_type",
                     "y0", "y1", "t1", "alpha", "rho") %in% names(result)))
  expect_true(all(result$antigen == "IpaB"))
})

test_that("prep_pop_params has correct number of rows", {
  model <- make_tiny_model()
  result <- prep_pop_params(model, "TestAg")
  # 4 iterations × 1 chain × 2 isotypes = 8
  expect_equal(nrow(result), 8)
})

test_that("prep_pop_params renames shape to rho", {
  model <- make_tiny_model()
  result <- prep_pop_params(model, "X")
  expect_true("rho" %in% names(result))
  expect_false("shape" %in% names(result))
})

# ── transform_pop_params ────────────────────────────────────────────────────

test_that("transform_pop_params adds natural-scale columns", {
  model <- make_tiny_model()
  draws <- prep_pop_params(model, "IpaB")
  result <- transform_pop_params(draws)

  new_cols <- c("y0_natural", "y1_natural", "t1_natural",
                "alpha_year", "rho_natural")
  expect_true(all(new_cols %in% names(result)))
})

test_that("transform_pop_params applies correct transformations", {
  df <- tibble::tibble(
    Iteration = 1L, Chain = 1L, antigen = "X", Iso_type = "IgG",
    y0 = log(100), y1 = log(900), t1 = log(5),
    alpha = log(0.01), rho = log(0.5)
  )
  result <- transform_pop_params(df)

  expect_equal(result$y0_natural, 100, tolerance = 1e-6)
  expect_equal(result$y1_natural, 100 + 900, tolerance = 1e-6)
  expect_equal(result$t1_natural, 5, tolerance = 1e-6)
  expect_equal(result$alpha_year, 0.01 * 365, tolerance = 1e-6)
  expect_equal(result$rho_natural, 0.5 + 1, tolerance = 1e-6)
})

# ── summarise_pop_params ────────────────────────────────────────────────────

test_that("summarise_pop_params returns one row per antigen-isotype", {
  model <- make_tiny_model()
  result <- prep_pop_params(model, "IpaB") |>
    transform_pop_params() |>
    summarise_pop_params()

  expect_equal(nrow(result), 2)  # IgG + IgA
  expect_true(all(c("y0_med", "y0_lo", "y0_hi") %in% names(result)))
})

# ── fmt_mci ─────────────────────────────────────────────────────────────────

test_that("fmt_mci formats correctly", {
  result <- fmt_mci(1.23, 0.45, 2.01, digits = 2)
  expect_type(result, "character")
  expect_match(result, "1.23")
  expect_match(result, "0.45")
})

test_that("fmt_mci supports scientific notation", {
  result <- fmt_mci(3.2e-9, 1.1e-10, 8.7e-8, digits = 2, sci = TRUE)
  expect_match(result, "e")
})

# ── format_param_table ──────────────────────────────────────────────────────

test_that("format_param_table produces Biomarker column", {
  model <- make_tiny_model()
  result <- prep_pop_params(model, "IpaB") |>
    transform_pop_params() |>
    summarise_pop_params() |>
    format_param_table()

  expect_true("Biomarker" %in% names(result))
  expect_true(all(grepl("IpaB", result$Biomarker)))
  # All param columns should be character (formatted strings)
  expect_type(result$y0, "character")
})
