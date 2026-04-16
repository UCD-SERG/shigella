#' Build the All-Models Parameter Comparison Table (S4 Table)
#'
#' Combines kinetic parameter summaries from the overall, serotype-specific,
#' and combined *S. flexneri* models for *S. flexneri* 2a, 3a, and *S. sonnei*,
#' and formats them into a publication-ready \pkg{gt} table.  Corresponds to
#' S4 Table in the Chapter 1 manuscript.
#'
#' @param pop_overall A summarised tibble (from [summarise_pop_params()])
#'   for the overall model (n = 48), containing all five antigens.
#' @param pop_serospec A summarised tibble for the serotype-specific models
#'   (Sf2a, Sf3a, Sonnei).
#' @param pop_combined A summarised tibble for the combined *S. flexneri*
#'   model (Sf2a, Sf3a), or `NULL` to omit.
#' @param n_by_model Optional named integer vector or data frame specifying
#'   sample sizes by antigen–model combination.  If `NULL`, hard-coded
#'   defaults from the SOSAR cohort are used.
#'
#' @return A `gt_tbl` object.
#'
#' @details
#' Only the antigen–isotype rows that appear in `pop_serospec` or
#' `pop_combined` are included, supplemented by the corresponding overall-model
#' rows.  IpaB and Sf6 (overall-only) are intentionally excluded to keep this
#' table focused on the model comparison.
#'
#' @examples
#' \dontrun{
#' gt_s4 <- build_table_s4_gt(
#'   pop_overall  = pop_table_sum,
#'   pop_serospec = pop_table_sum2,
#'   pop_combined = pop_table_sum3
#' )
#' print(gt_s4)
#' }
#'
#' @importFrom dplyr bind_rows filter mutate select case_when arrange transmute
#' @importFrom gt gt fmt_number cols_label cols_width tab_header tab_source_note
#'   tab_options md px
#' @export
build_table_s4_gt <- function(pop_overall,
                               pop_serospec,
                               pop_combined = NULL,
                               n_by_model   = NULL) {

  # Default sample-size lookup (SOSAR cohort)
  if (is.null(n_by_model)) {
    n_by_model <- data.frame(
      antigen = c("Sf2a", "Sf2a", "Sf2a",
                  "Sonnei", "Sonnei",
                  "Sf3a", "Sf3a", "Sf3a"),
      model   = c("Overall", "Serotype-specific", "Combined flexneri",
                  "Overall", "Serotype-specific",
                  "Overall", "Serotype-specific", "Combined flexneri"),
      n       = c(48L, 17L, 25L,
                  48L, 11L,
                  48L,  8L, 25L),
      stringsAsFactors = FALSE
    )
  }

  # Antigens for which a comparison exists
  compare_antigens <- c("Sf2a", "Sf3a", "Sonnei")

  overall_for_supp <- pop_overall |>
    dplyr::filter(.data$antigen %in% compare_antigens) |>
    dplyr::mutate(model = "Overall")

  sero_for_supp <- pop_serospec |>
    dplyr::mutate(model = "Serotype-specific")

  all_rows <- dplyr::bind_rows(
    overall_for_supp,
    sero_for_supp,
    if (!is.null(pop_combined))
      dplyr::mutate(pop_combined, model = "Combined flexneri")
  ) |>
    dplyr::mutate(
      antigen  = factor(.data$antigen,
                        levels = c("Sf2a", "Sonnei", "Sf3a")),
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA")),
      model    = factor(.data$model,
                        levels = c("Overall",
                                   "Serotype-specific",
                                   "Combined flexneri"))
    ) |>
    dplyr::arrange(.data$antigen, .data$Iso_type, .data$model)

  # Attach sample sizes
  all_rows <- dplyr::left_join(
    all_rows,
    n_by_model,
    by = c("antigen", "model")
  )

  fmt_f <- function(med, lo, hi, digits = 1) {
    sprintf(
      "%s (%s\u2013%s)",
      formatC(med, format = "f", digits = digits),
      formatC(lo,  format = "f", digits = digits),
      formatC(hi,  format = "f", digits = digits)
    )
  }

  display_df <- all_rows |>
    dplyr::transmute(
      Antigen          = .data$antigen,
      Isotype          = .data$Iso_type,
      Model            = .data$model,
      n                = .data$n,
      `y0 (baseline)`  = fmt_f(.data$y0_med,    .data$y0_lo,    .data$y0_hi),
      `y1 (peak)`      = fmt_f(.data$y1_med,    .data$y1_lo,    .data$y1_hi, 0),
      `t1 (days)`      = fmt_f(.data$t1_med,    .data$t1_lo,    .data$t1_hi),
      `alpha (decay)`  = fmt_f(.data$alpha_med, .data$alpha_lo, .data$alpha_hi, 2),
      `rho (shape)`    = fmt_f(.data$rho_med,   .data$rho_lo,   .data$rho_hi, 2)
    )

  display_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md(
        "**Population-level kinetic parameter estimates across all modelling approaches**"
      )
    ) |>
    gt::cols_label(
      Antigen         = "Antigen",
      Isotype         = "Isotype",
      Model           = "Model",
      n               = "n",
      `y0 (baseline)` = gt::md("y<sub>0</sub> (baseline)"),
      `y1 (peak)`     = gt::md("y<sub>1</sub> (peak)"),
      `t1 (days)`     = gt::md("t<sub>1</sub> (days)"),
      `alpha (decay)` = gt::md("&alpha; (decay)"),
      `rho (shape)`   = gt::md("&rho; (shape)")
    ) |>
    gt::tab_source_note(gt::md(
      "Values are posterior medians (95% credible intervals)."
    )) |>
    gt::tab_source_note(gt::md(
      "The combined flexneri model pools *S. flexneri* 2a- and 3a-infected individuals (*n* = 25)."
    )) |>
    gt::tab_source_note(gt::md(
      "*S. sonnei* was modelled as overall and serotype-specific only."
    )) |>
    gt::tab_source_note(gt::md(
      "MCMC: 4 chains, 25,000 adaptation, 50,000 burn-in, 15,000 retained samples per chain."
    )) |>
    gt::cols_width(
      Antigen         ~ gt::px(90),
      Isotype         ~ gt::px(80),
      Model           ~ gt::px(150),
      n               ~ gt::px(50),
      gt::everything() ~ gt::px(130)
    ) |>
    gt::tab_options(
      table.font.size        = gt::px(12),
      source_notes.font.size = gt::px(10)
    )
}


#' Build the Raw Per-Individual MAE Table (S5 Table)
#'
#' Creates a wide-format tibble (or kable-formatted HTML table) of
#' per-individual MAE on the log-MFI scale, joining overall, serotype-specific,
#' and combined-flexneri MAE results.  Corresponds to S5 Table in the
#' Chapter 1 manuscript.
#'
#' @param mae_overall A tibble from [get_mae()] for the overall model.
#' @param mae_serospec A tibble from [get_mae()] for serotype-specific models,
#'   or `NULL`.
#' @param mae_combined A tibble from [get_mae()] for the combined flexneri
#'   model, or `NULL`.
#' @param as_kable Logical.  If `TRUE` (default `FALSE`), return a
#'   `kableExtra`-formatted HTML table with scrollable box.  If `FALSE`,
#'   return a plain tibble.
#'
#' @return Either a tibble or a `kableExtra` HTML table, depending on
#'   `as_kable`.
#'
#' @examples
#' \dontrun{
#' tbl <- build_table_s5(
#'   mae_overall  = mae_overall,
#'   mae_serospec = mae_serospec,
#'   mae_combined = mae_combined,
#'   as_kable     = TRUE
#' )
#' }
#'
#' @importFrom dplyr full_join rename mutate arrange across starts_with
#' @importFrom tibble tibble
#' @export
build_table_s5 <- function(mae_overall,
                            mae_serospec = NULL,
                            mae_combined = NULL,
                            as_kable     = FALSE) {

  out <- mae_overall |>
    dplyr::rename(MAE_overall = "mae")

  if (!is.null(mae_serospec)) {
    out <- dplyr::full_join(
      out,
      mae_serospec |> dplyr::rename(MAE_serospec = "mae"),
      by = c("sid", "antigen", "Iso_type")
    )
  } else {
    out$MAE_serospec <- NA_real_
  }

  if (!is.null(mae_combined)) {
    out <- dplyr::full_join(
      out,
      mae_combined |> dplyr::rename(MAE_combined = "mae"),
      by = c("sid", "antigen", "Iso_type")
    )
  } else {
    out$MAE_combined <- NA_real_
  }

  out <- out |>
    dplyr::mutate(
      best = dplyr::case_when(
        is.na(.data$MAE_serospec) & is.na(.data$MAE_combined) ~ "Overall only",
        is.na(.data$MAE_serospec)                              ~ "Overall only",
        .data$MAE_serospec < .data$MAE_overall                 ~ "Sero-specific",
        .data$MAE_serospec > .data$MAE_overall                 ~ "Overall",
        TRUE                                                   ~ "Tie"
      )
    ) |>
    dplyr::arrange(.data$antigen, .data$Iso_type, .data$sid) |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("MAE"), ~ round(.x, 4))
    )

  if (!as_kable) return(out)

  # kableExtra HTML table (requires knitr + kableExtra)
  if (!requireNamespace("knitr",      quietly = TRUE) ||
      !requireNamespace("kableExtra", quietly = TRUE)) {
    warning("'knitr' and 'kableExtra' are required for as_kable = TRUE.")
    return(out)
  }

  out |>
    knitr::kable(
      caption   = "Per-individual MAE on the log-MFI scale for all models.",
      col.names = c(
        "ID", "Antigen", "Isotype",
        "MAE (Overall)", "MAE (Sero-specific)", "MAE (Combined)", "Better"
      )
    ) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "condensed"),
      font_size         = 9,
      full_width        = TRUE
    ) |>
    kableExtra::scroll_box(height = "500px")
}
