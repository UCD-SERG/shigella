# tests/testthat/test-figures.R
# Tests for build_figure4_row, fig5_with_individuals, fig5_ipab_with_age

# ── Helpers ─────────────────────────────────────────────────────────────────

make_mini_model <- function(sids, n_iter = 5L) {
  params <- c("y0", "y1", "t1", "alpha", "shape")
  base <- c(5.0, 10.0, -0.5, -15.0, 0.3)

  indiv_rows <- list()
  pop_rows   <- list()

  for (sid in sids) {
    for (iso in c("IgG", "IgA")) {
      for (it in seq_len(n_iter)) {
        vals <- base + stats::rnorm(5, 0, 0.1)
        for (p in seq_along(params)) {
          indiv_rows[[length(indiv_rows) + 1L]] <- tibble::tibble(
            Iteration = it, Chain = 1L, Parameter = params[p],
            Iso_type = iso, Stratification = "None",
            Subject = sid, value = vals[p]
          )
        }
      }
    }
  }

  for (iso in c("IgG", "IgA")) {
    for (it in seq_len(n_iter)) {
      vals <- base + stats::rnorm(5, 0, 0.05)
      for (p in seq_along(params)) {
        pop_rows[[length(pop_rows) + 1L]] <- tibble::tibble(
          Iteration = it, Chain = 1L, Parameter = params[p],
          Iso_type = iso, Stratification = "None",
          Population_Parameter = "mu.par", value = vals[p]
        )
      }
      # prec.logy
      pop_rows[[length(pop_rows) + 1L]] <- tibble::tibble(
        Iteration = it, Chain = 1L, Parameter = "prec.logy",
        Iso_type = iso, Stratification = "None",
        Population_Parameter = "prec.logy", value = stats::runif(1, 5, 15)
      )
    }
  }

  structure(
    dplyr::bind_rows(indiv_rows),
    class = c("sr_model", "tbl_df", "tbl", "data.frame"),
    population_params = dplyr::bind_rows(pop_rows)
  )
}

make_mini_casedata <- function(sids) {
  rows <- list()
  for (sid in sids) {
    for (iso in c("IgG", "IgA")) {
      for (tp in c(2, 30, 90)) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          id = sid, antigen_iso = iso, timeindays = tp,
          result = exp(stats::runif(1, 6, 11)),
          sid = sid, isotype_name = iso, cohort_name = "Sf2a"
        )
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  structure(
    out,
    class = c("case_data", "tbl_df", "tbl", "data.frame"),
    id_var = "id", biomarker_var = "antigen_iso",
    timeindays = "timeindays", value_var = "result"
  )
}

# ── build_figure4_row ───────────────────────────────────────────────────────

test_that("build_figure4_row returns a patchwork object", {
  sids <- c("S1", "S2")
  model_ov   <- make_mini_model(sids)
  model_sero <- make_mini_model("S1")
  data_ov    <- make_mini_casedata(sids)
  data_sero  <- make_mini_casedata("S1")

  result <- build_figure4_row(
    sid = "S1",
    model_overall = model_ov,
    model_sero    = model_sero,
    data_overall  = data_ov,
    data_sero     = data_sero,
    row_label     = "Test row",
    antigen_label = "TestAg",
    times         = seq(0, 90, by = 30)
  )
  # patchwork objects inherit from gg
  expect_true(inherits(result, "patchwork") || inherits(result, "gg"))
})

# ── fig5_with_individuals ───────────────────────────────────────────────────

test_that("fig5_with_individuals returns a ggplot", {
  model <- make_mini_model(c("S1", "S2"))
  cd    <- make_mini_casedata(c("S1", "S2"))

  result <- fig5_with_individuals(
    model_obj     = model,
    raw_dataset   = cd,
    antigen_label = "Test antigen",
    t_grid        = seq(0, 90, by = 30)
  )
  expect_s3_class(result, "gg")
})

# ── fig5_ipab_with_age ──────────────────────────────────────────────────────

test_that("fig5_ipab_with_age returns a ggplot", {
  model_all   <- make_mini_model(c("S1", "S2", "S3"))
  model_u5    <- make_mini_model(c("S1"))
  model_p5    <- make_mini_model(c("S2", "S3"))
  cd_all      <- make_mini_casedata(c("S1", "S2", "S3"))
  cd_u5       <- make_mini_casedata(c("S1"))
  cd_p5       <- make_mini_casedata(c("S2", "S3"))

  result <- fig5_ipab_with_age(
    model_overall = model_all,
    model_under5  = model_u5,
    model_plus5   = model_p5,
    raw_overall   = cd_all,
    raw_under5    = cd_u5,
    raw_plus5     = cd_p5,
    t_grid        = seq(0, 90, by = 30)
  )
  expect_s3_class(result, "gg")
})

# ── .extract_individual_trajectories (internal helper) ──────────────────────

test_that(".extract_individual_trajectories returns correct columns", {
  cd <- make_mini_casedata(c("S1", "S2"))
  result <- shigella:::.extract_individual_trajectories(cd)

  expect_true(all(c("id", "t", "value", "iso") %in% names(result)))
  expect_s3_class(result$iso, "factor")
  expect_equal(levels(result$iso), c("IgG", "IgA"))
})
