#' S1 Table: prior distributions for the population-level parameters
#'
#' The hierarchical two-phase model places priors on the log-scale population
#' means, the between-individual precision (Wishart), and the observation
#' precision (Gamma). Values are fixed (the primary specification) and match
#' `prior_settings` in `data-raw/_config.R`.
#'
#' @return A `gt_tbl`.
#' @export
table_s1_priors <- function() {
  table_s1 <- tibble::tribble(
    ~Level,            ~Parameter,       ~Distribution,        ~Hyperparameters,               ~Interpretation, # nolint: line_length_linter.
    "Population mean", "\u03BC(log y\u2080)",     "Normal(m, \u03C4\u207B\u00B9)", "m = 8, \u03C4 = 0.25",         "Baseline MFI ~3000 (exp(8))", # nolint: line_length_linter.
    "Population mean", "\u03BC(log y\u2081)",     "Normal(m, \u03C4\u207B\u00B9)", "m = 10, \u03C4 = 1\u00D710\u207B\u2074", "Peak MFI ~22000 (exp(10))", # nolint: line_length_linter.
    "Population mean", "\u03BC(log t\u2081)",     "Normal(m, \u03C4\u207B\u00B9)", "m = 0.5, \u03C4 = 0.5",        "Time to peak ~1.6 days (exp(0.5))", # nolint: line_length_linter.
    "Population mean", "\u03BC(log \u03B1)",      "Normal(m, \u03C4\u207B\u00B9)", "m = -4, \u03C4 = 0.001",       "Decay rate ~0.018 yr\u207B\u00B9 (exp(-4))", # nolint: line_length_linter.
    "Population mean", "\u03BC(log(\u03C1-1))",   "Normal(m, \u03C4\u207B\u00B9)", "m = -1, \u03C4 = 0.25",        "Decay shape ~1.37 (1+exp(-1))", # nolint: line_length_linter.
    "Covariance",      "\u03A3\u2C7C\u207B\u00B9", "Wishart(\u03A9, \u03BD)",     "\u03A9 = diag(1,1,1,1,1), \u03BD = 8", "5\u00D75 precision; weakly informative", # nolint: line_length_linter.
    "Measurement",     "\u03C4\u2C7C (precision)", "Gamma(a, b)",                "a = 1, b = 1",                 "Observation noise; weakly informative" # nolint: line_length_linter.
  )

  table_s1 |>
    gt::gt(rowname_col = "Parameter", groupname_col = "Level") |>
    gt::tab_header(
      title = gt::md("**Prior distributions for population-level parameters**")
    ) |>
    gt::cols_label(
      Distribution    = "Distribution",
      Hyperparameters = "Hyperparameters",
      Interpretation  = "Interpretation"
    ) |>
    gt::tab_row_group(
      label = "Population-level means (\u03BC\u2C7C)",
      rows  = .data$Level == "Population mean"
    ) |>
    gt::tab_row_group(
      label = "Between-individual covariance",
      rows  = .data$Level == "Covariance"
    ) |>
    gt::tab_row_group(
      label = "Observation model",
      rows  = .data$Level == "Measurement"
    ) |>
    gt::cols_width(
      Parameter ~ gt::px(140),
      Distribution ~ gt::px(140),
      Hyperparameters ~ gt::px(180),
      Interpretation ~ gt::px(240)
    ) |>
    gt::tab_source_note(
      source_note = gt::md("All priors were specified on the logarithmic scale. Precision (&tau;) is the inverse of variance.") # nolint: line_length_linter.
    ) |>
    gt::tab_source_note(
      source_note = gt::md("Hyperparameter values were informed by prior predictive simulations and existing applications of the two-phase model to enteric pathogens.") # nolint: line_length_linter.
    ) |>
    gt::tab_options(
      table.font.size = gt::px(11),
      source_notes.font.size = gt::px(10)
    )
}
