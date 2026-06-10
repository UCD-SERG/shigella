# ── Internal helpers ────────────────────────────────

# Wide display frame for S3: one row per biomarker x prior, with delta-MAE and a
# Serotype/Overall/Similar conclusion + a markdown group label.
#' @keywords internal
#' @noRd
.sensitivity_wide <- function(sensitivity_results) {
  sensitivity_results |>
    tidyr::pivot_wider(names_from = "model", values_from = "mae",
                       names_prefix = "mae_") |>
    dplyr::mutate(
      delta_mae = .data$mae_serotype - .data$mae_overall,
      conclusion = dplyr::case_when(
        .data$delta_mae < -0.1 ~ "Serotype better",
        .data$delta_mae >  0.1 ~ "Overall better",
        TRUE                   ~ "Similar"
      ),
      biomarker_group = dplyr::case_when(
        .data$biomarker == "Sonnei IgG" ~ "*S. sonnei* IgG",
        .data$biomarker == "Sonnei IgA" ~ "*S. sonnei* IgA",
        .data$biomarker == "Sf3a IgG"   ~ "*S. flexneri* 3a IgG",
        .data$biomarker == "Sf3a IgA"   ~ "*S. flexneri* 3a IgA",
        TRUE ~ as.character(.data$biomarker)
      )
    ) |>
    dplyr::select(
      "biomarker_group",
      Biomarker = "biomarker", Prior = "prior",
      `MAE (Overall)` = "mae_overall", `MAE (Serotype)` = "mae_serotype",
      `Delta MAE` = "delta_mae", Conclusion = "conclusion"
    )
}

# ── Main function ─────────────────────────────────

#' S3 Table: prior-specification robustness (overall vs serotype MAE)
#'
#' Data shaping is in [.sensitivity_wide()]; this adds gt styling.
#'
#' @param sensitivity_results Output of [build_sensitivity_results()].
#' @return A `gt_tbl`.
#' @export
table_s3_sensitivity <- function(sensitivity_results) {
  .sensitivity_wide(sensitivity_results) |>
    gt::gt(groupname_col = "biomarker_group") |>
    gt::fmt_number(columns = c("MAE (Overall)", "MAE (Serotype)", "Delta MAE"),
                   decimals = 2) |>
    gt::cols_label(
      Biomarker = "Biomarker", Prior = "Prior",
      `MAE (Overall)` = "MAE (Overall)", `MAE (Serotype)` = "MAE (Serotype)",
      `Delta MAE` = gt::md("&Delta;MAE"), Conclusion = "Conclusion"
    ) |>
    gt::cols_width(
      Biomarker ~ gt::px(120), Prior ~ gt::px(100),
      `MAE (Overall)` ~ gt::px(110), `MAE (Serotype)` ~ gt::px(115),
      `Delta MAE` ~ gt::px(80), Conclusion ~ gt::px(130)
    ) |>
    gt::tab_header(
      title = gt::md("**Sensitivity analysis: prior specification robustness**")
    ) |>
    gt::tab_source_note(gt::md(paste0(
      "For each biomarker, models were refit under three prior configurations: ",
      "primary (current), diffuse (2\u00D7 SD), and informative (0.5\u00D7 SD)."))) |>
    gt::tab_source_note(gt::md(paste0(
      "&Delta;MAE = MAE(serotype) &minus; MAE(overall); negative values indicate ",
      "serotype-specific advantage."))) |>
    gt::tab_source_note(gt::md(paste0(
      "The pattern&mdash;*S. sonnei* favoring serotype-specific and ",
      "*S. flexneri* 3a favoring overall&mdash;is consistent across all prior ",
      "configurations, indicating that model comparison conclusions are driven ",
      "by data structure and cross-reactivity, not prior choice."))) |>
    gt::tab_options(
      table.font.size = gt::px(10),
      source_notes.font.size = gt::px(10)
    )
}
