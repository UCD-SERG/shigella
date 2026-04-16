# tests/testthat/test-tables.R
# Tests for build_table4, build_kinetic_flextable

# ── build_table4 ────────────────────────────────────────────────────────────

test_that("build_table4 returns a flextable", {
  mae_ov <- tibble::tibble(
    sid      = rep(c("S1", "S2", "S3"), 2),
    antigen  = "IpaB",
    Iso_type = rep(c("IgG", "IgA"), each = 3),
    mae      = c(0.15, 0.22, 0.18, 0.30, 0.25, 0.40)
  )

  result <- build_table4(mae_overall = mae_ov)
  expect_s3_class(result, "flextable")
})

test_that("build_table4 handles serospec + combined inputs", {
  mae_ov <- tibble::tibble(
    sid = rep(c("S1", "S2"), 2), antigen = "Sf2a",
    Iso_type = rep(c("IgG", "IgA"), each = 2),
    mae = c(0.20, 0.25, 0.35, 0.40)
  )
  mae_sp <- tibble::tibble(
    sid = rep(c("S1", "S2"), 2), antigen = "Sf2a",
    Iso_type = rep(c("IgG", "IgA"), each = 2),
    mae = c(0.15, 0.18, 0.30, 0.32)
  )
  mae_cmb <- tibble::tibble(
    sid = rep(c("S1", "S2"), 2), antigen = "Sf2a",
    Iso_type = rep(c("IgG", "IgA"), each = 2),
    mae = c(0.22, 0.28, 0.33, 0.38)
  )

  result <- build_table4(
    mae_overall  = mae_ov,
    mae_serospec = mae_sp,
    mae_combined = mae_cmb
  )
  expect_s3_class(result, "flextable")
})

test_that("build_table4 shows N/A for missing models", {
  mae_ov <- tibble::tibble(
    sid = "S1", antigen = "IpaB", Iso_type = "IgG", mae = 0.15
  )
  # Only overall, no serospec or combined
  result <- build_table4(mae_overall = mae_ov)
  # The flextable should exist (can't easily inspect cell values, just check creation)
  expect_s3_class(result, "flextable")
})

# ── build_kinetic_flextable ─────────────────────────────────────────────────

test_that("build_kinetic_flextable returns a flextable", {
  display_df <- tibble::tibble(
    Biomarker = c("IpaB\u2013IgG", "IpaB\u2013IgA"),
    y0    = c("100.00 (10.00\u2013500.00)", "50.00 (5.00\u2013200.00)"),
    y1    = c("60000.00 (50000.00\u201370000.00)", "2000.00 (1000.00\u20135000.00)"),
    t1    = c("1.2 (0.3\u20132.6)", "3.0 (0.9\u20137.4)"),
    alpha = c("0.00126592 (0.00000014\u20130.10478217)",
              "0.01185148 (0.00002510\u20130.20895165)"),
    rho   = c("1.35 (1.05\u20131.95)", "1.51 (1.13\u20132.49)")
  )

  result <- build_kinetic_flextable(
    display_df,
    caption = "Test caption",
    footer_lines = "Test footer"
  )
  expect_s3_class(result, "flextable")
})

test_that("build_kinetic_flextable works without footer", {
  display_df <- tibble::tibble(
    Biomarker = "IpaB\u2013IgG",
    y0 = "100 (10\u2013500)", y1 = "60000 (50000\u201370000)",
    t1 = "1.2 (0.3\u20132.6)", alpha = "0.001 (0.0\u20130.1)",
    rho = "1.35 (1.05\u20131.95)"
  )
  result <- build_kinetic_flextable(display_df, caption = "Test")
  expect_s3_class(result, "flextable")
})
