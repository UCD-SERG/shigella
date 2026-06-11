# ── Internal helpers ────────────────────────────────

# Format a kinetic-parameter summary (one row per Biomarker) for display.
#' @keywords internal
#' @noRd
.format_kinetic_rows <- function(summary_tbl, extra_cols = character(0)) {
  summary_tbl |>
    dplyr::mutate(
      Biomarker = paste0(as.character(.data$antigen), "\u2013", as.character(.data$Iso_type)),
      y0    = fmt_mci(.data$y0_med, .data$y0_lo, .data$y0_hi, digits = 2),
      y1    = fmt_mci(.data$y1_med, .data$y1_lo, .data$y1_hi, digits = 2),
      t1    = fmt_mci(.data$t1_med, .data$t1_lo, .data$t1_hi, digits = 1),
      alpha = fmt_mci(.data$alpha_med, .data$alpha_lo, .data$alpha_hi, digits = 8),
      rho   = fmt_mci(.data$rho_med, .data$rho_lo, .data$rho_hi, digits = 2)
    ) |>
    dplyr::select(dplyr::all_of(c("Biomarker", extra_cols,
                                  "y0", "y1", "t1", "alpha", "rho")))
}

# Pull the parameter-summary row for one biomarker from whichever model won,
# adding a display label with that model's n.
#' @keywords internal
#' @noRd
.params_for_best_row <- function(antigen_val, iso_val, model_val,
                                 sum_overall, sum_serospec, sum_combined) {
  if (model_val == "Overall") {
    src <- sum_overall; model_display <- "Overall (n=48)"
  } else if (model_val == "Serotype-specific") {
    src <- sum_serospec
    n_val <- dplyr::case_when(antigen_val == "Sf2a" ~ 17L, antigen_val == "Sf3a" ~ 8L,
                              antigen_val == "Sonnei" ~ 11L, TRUE ~ NA_integer_)
    model_display <- sprintf("Serotype-specific (n=%d)", n_val)
  } else {
    src <- sum_combined; model_display <- "Combined flexneri (n=25)"
  }
  src |>
    dplyr::filter(.data$antigen == antigen_val, .data$Iso_type == iso_val) |>
    dplyr::mutate(Model = model_display)
}

# flextable styling + footnotes for Table 2.
#' @keywords internal
#' @noRd
.table2_flextable <- function(best_print) {
  best_print |>
    flextable::flextable() |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 7, part = "body") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1:2, align = "left", part = "body") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption(paste0(
      "Antibody kinetic parameters from the best-performing model for each ",
      "antigen-isotype combination.")) |>
    flextable::add_footer_lines(values = c(
      paste0("Values are posterior medians (95% credible intervals) for the ",
             "population mean (mu.par) of each parameter under the best-performing model."),
      paste0("Best model was determined by per-individual MAE comparisons on shared ",
             "individuals. Parameters: y0, baseline; y1, peak; t1, time to peak; ",
             "alpha, decay rate (year^-1); rho, decay shape."))) |>
    flextable::fontsize(size = 6, part = "footer") |>
    flextable::italic(part = "footer") |>
    flextable::color(color = "gray40", part = "footer")
}

# ── Main function ─────────────────────────────────

#' Table 2: kinetic parameters from the best-performing model per biomarker
#'
#' For each biomarker the winning model's parameter row is pulled by
#' `.params_for_best_row()`, formatted, then styled by `.table2_flextable()`.
#'
#' @param best_lookup Output of [select_best_models()].
#' @param sum_overall,sum_serospec,sum_combined Parameter summaries
#'   ([pop_param_summary()] output) for each model class.
#' @return A flextable.
#' @export
table2_kinetic_params <- function(best_lookup, sum_overall, sum_serospec, sum_combined) {
  best_data <- purrr::pmap_dfr(best_lookup, function(antigen, Iso_type, best_model)
    .params_for_best_row(antigen, Iso_type, best_model,
                         sum_overall, sum_serospec, sum_combined))

  best_print <- .format_kinetic_rows(best_data, extra_cols = "Model")
  .add_kinetic_headers(.table2_flextable(best_print)) # nolint: object_usage_linter
}
