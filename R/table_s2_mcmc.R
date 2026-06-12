#' S2 Table: MCMC settings
#'
#' @return A `gt_tbl`.
#' @export
table_s2_mcmc <- function() {
  table_s2 <- tibble::tribble(
    ~Setting,                     ~Value,       ~Description,
    "Number of chains",           "4",          "Independent chains with dispersed initial values", # nolint: line_length_linter.
    "Adaptation iterations",      "25,000",     "JAGS sampler tuning period",
    "Burn-in iterations",         "50,000",     "Discarded to ensure convergence", # nolint: line_length_linter.
    "Retained samples per chain", "15,000",     "Post-thinning samples",
    "Total posterior samples",    "60,000",     "4 chains \u00D7 15,000 = 60,000", # nolint: line_length_linter.
    "Software",                   "JAGS 4.3.4", "Via rjags/runjags in R 4.4.3",
    "Convergence criterion",      "R\u0302 < 1.1", "Gelman-Rubin diagnostic for all monitored parameters" # nolint: line_length_linter.
  )

  table_s2 |>
    gt::gt() |>
    gt::tab_header(title = gt::md("**MCMC settings**")) |>
    gt::cols_label(Setting = "Setting", Value = "Value", Description = "Description") |> # nolint: line_length_linter.
    gt::cols_width(
      Setting ~ gt::px(180),
      Value ~ gt::px(110),
      Description ~ gt::px(260)
    ) |>
    gt::tab_source_note(
      source_note = gt::md("Markov chain Monte Carlo (MCMC) settings for all models.") # nolint: line_length_linter.
    ) |>
    gt::tab_source_note(
      source_note = gt::md("Convergence was assessed using visual inspection of trace plots and the potential scale reduction factor (R&#770;). All monitored parameters achieved R&#770; < 1.05.") # nolint: line_length_linter.
    ) |>
    gt::tab_options(
      table.font.size = gt::px(11),
      source_notes.font.size = gt::px(10)
    )
}
