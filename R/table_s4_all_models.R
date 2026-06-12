# Internal helpers

# Stack the three model summaries (overall Sf2a/Sf3a/Sonnei rows + serotype +
# combined flexneri) into one factor-ordered frame.
#' @keywords internal
#' @noRd
.combine_model_summaries <- function(sum_overall, sum_serospec, sum_combined) {
  overall_for_supp <- sum_overall |>
    dplyr::filter(.data$antigen %in% c("Sf2a", "Sf3a", "Sonnei")) |>
    dplyr::mutate(model = "Overall")
  sero_for_supp <- sum_serospec |> dplyr::mutate(model = "Serotype-specific")
  comb_for_supp <- sum_combined |> dplyr::mutate(model = "Combined flexneri")

  dplyr::bind_rows(overall_for_supp, sero_for_supp, comb_for_supp) |>
    dplyr::mutate(
      antigen  = factor(.data$antigen, levels = c("Sf2a", "Sonnei", "Sf3a")),
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA")),
      model    = factor(.data$model,
                        levels = c("Overall", "Serotype-specific", "Combined flexneri")) # nolint: line_length_linter.
    ) |>
    dplyr::arrange(.data$antigen, .data$Iso_type, .data$model)
}

# Attach the per-cell n and format each parameter as "median (lo, hi)".
#' @keywords internal
#' @noRd
.s4_display <- function(supp_all) {
  supp_all |>
    dplyr::mutate(
      n = dplyr::case_when(
        .data$antigen == "Sf2a"   & .data$model == "Overall"            ~ 48L,
        .data$antigen == "Sf2a"   & .data$model == "Serotype-specific"  ~ 17L,
        .data$antigen == "Sf2a"   & .data$model == "Combined flexneri"  ~ 25L,
        .data$antigen == "Sonnei" & .data$model == "Overall"            ~ 48L,
        .data$antigen == "Sonnei" & .data$model == "Serotype-specific"  ~ 11L,
        .data$antigen == "Sf3a"   & .data$model == "Overall"            ~ 48L,
        .data$antigen == "Sf3a"   & .data$model == "Serotype-specific"  ~ 8L,
        .data$antigen == "Sf3a"   & .data$model == "Combined flexneri"  ~ 25L,
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::transmute(
      Antigen = .data$antigen, Isotype = .data$Iso_type, Model = .data$model, n = .data$n, # nolint: line_length_linter.
      `y0 (baseline)` = fmt_mci(.data$y0_med, .data$y0_lo, .data$y0_hi, digits = 1), # nolint: line_length_linter.
      `y1 (peak)`     = fmt_mci(.data$y1_med, .data$y1_lo, .data$y1_hi, digits = 0), # nolint: line_length_linter.
      `t1 (days)`     = fmt_mci(.data$t1_med, .data$t1_lo, .data$t1_hi, digits = 1), # nolint: line_length_linter.
      `alpha (decay)` = fmt_mci(.data$alpha_med, .data$alpha_lo, .data$alpha_hi, digits = 2), # nolint: line_length_linter.
      `rho (shape)`   = fmt_mci(.data$rho_med, .data$rho_lo, .data$rho_hi, digits = 2) # nolint: line_length_linter.
    )
}

# Main function

#' S4 Table: population kinetic parameters across all modeling approaches
#'
#' Composes `.combine_model_summaries()` and `.s4_display()`, then gt styling.
#' The three inputs are the `pop_param_summary()` outputs from the main
#' manuscript (`pop_table_sum`, `pop_table_sum2`, `pop_table_sum3`).
#'
#' @param sum_overall `pop_param_summary()` table for the overall model
#'   (all five antigens; only Sf2a/Sf3a/Sonnei rows are used here).
#' @param sum_serospec `pop_param_summary()` table for the serotype-specific
#'   models (Sf2a n=17, Sf3a n=8, Sonnei n=11).
#' @param sum_combined `pop_param_summary()` table for the combined flexneri
#'   model (Sf2a/Sf3a, n=25).
#' @return A `gt_tbl`.
#' @export
table_s4_all_models <- function(sum_overall, sum_serospec, sum_combined) {
  .combine_model_summaries(sum_overall, sum_serospec, sum_combined) |>
    .s4_display() |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Population-level kinetic parameter estimates across all modeling approaches**") # nolint: line_length_linter.
    ) |>
    gt::cols_label(
      Antigen = "Antigen", Isotype = "Isotype", Model = "Model", n = "n",
      `y0 (baseline)` = gt::md("y<sub>0</sub> (baseline)"),
      `y1 (peak)`     = gt::md("y<sub>1</sub> (peak)"),
      `t1 (days)`     = gt::md("t<sub>1</sub> (days)"),
      `alpha (decay)` = gt::md("&alpha; (decay)"),
      `rho (shape)`   = gt::md("&rho; (shape)")
    ) |>
    gt::tab_source_note(gt::md("Values are posterior medians (95% credible intervals).")) |> # nolint: line_length_linter.
    gt::tab_source_note(gt::md("The combined flexneri model pools *S. flexneri* 2a- and 3a-infected individuals (*n* = 25).")) |> # nolint: line_length_linter.
    gt::tab_source_note(gt::md("*S. sonnei* was modeled as overall and serotype-specific only.")) |> # nolint: line_length_linter.
    gt::tab_source_note(gt::md("MCMC settings: 4 chains, 10,000,000 iterations, 50,000 burn-in.")) |> # nolint: line_length_linter.
    gt::cols_width(
      Antigen ~ gt::px(90), Isotype ~ gt::px(80), Model ~ gt::px(150),
      n ~ gt::px(50), gt::everything() ~ gt::px(130)
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      source_notes.font.size = gt::px(10)
    )
}
