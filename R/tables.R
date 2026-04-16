#' Build MAE Comparison Table Across Modeling Approaches
#'
#' Compares per-individual prediction accuracy (median absolute error, MAE)
#' across overall, serotype-specific, and combined *S. flexneri* models for
#' each antigen–isotype combination.  The best model is determined by counting
#' individual-level "wins" among participants shared across all compared models.
#'
#' @param mae_overall A tibble of per-individual MAE from the overall model,
#'   as returned by [get_mae()], with columns `sid`, `antigen`, `Iso_type`,
#'   `mae`.
#' @param mae_serospec A tibble of per-individual MAE from serotype-specific
#'   models, or `NULL`.
#' @param mae_combined A tibble of per-individual MAE from the combined
#'   *S. flexneri* model, or `NULL`.
#' @param n_lookup A tibble with columns `antigen` and `n` specifying the
#'   sample sizes to display.  If `NULL`, sample sizes are computed from
#'   the union of all input MAE tables.
#'
#' @return A `flextable` object ready for embedding in a Word document or
#'   Quarto PDF, with columns `Antigen`, `Isotype`, `n`, `Overall`,
#'   `Sero-specific`, `Combined flexneri`, `Best model`.  Values are
#'   `"median (Q1\u2013Q3)"` strings on the log-MFI scale.
#'
#' @details
#' The best model for each antigen–isotype is determined as follows:
#' \enumerate{
#'   \item Identify participants shared across all models being compared.
#'   \item For each shared participant, record which model yields the lowest
#'     MAE.
#'   \item The model with the most individual-level wins is declared the best.
#'   \item The `Best model` column reports the winning model name and
#'     `"wins / n_shared"` in parentheses.
#' }
#' When only one model is available for an antigen–isotype combination,
#' it is returned as the best model with the total number of participants.
#'
#' @examples
#' \dontrun{
#' tbl <- build_table4(
#'   mae_overall  = mae_overall,
#'   mae_serospec = mae_serospec,
#'   mae_combined = mae_combined,
#'   n_lookup     = n_lookup
#' )
#' tbl  # prints in RStudio viewer
#' }
#'
#' @importFrom dplyr bind_rows group_by summarise mutate left_join arrange
#'   filter select transmute distinct count slice_min ungroup n_distinct
#' @importFrom tidyr pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom flextable flextable theme_booktabs fontsize bold align
#'   set_table_properties set_caption
#' @export
build_table4 <- function(mae_overall,
                         mae_serospec = NULL,
                         mae_combined = NULL,
                         n_lookup     = NULL) {

  required_cols <- c("sid", "antigen", "Iso_type", "mae")

  standardize_mae <- function(df, model_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- tibble::as_tibble(df)
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(paste0("Input '", model_name, "' missing columns: ",
                  paste(missing_cols, collapse = ", ")))
    }
    df |>
      dplyr::transmute(
        sid      = as.character(.data$sid),
        antigen  = as.character(.data$antigen),
        Iso_type = as.character(.data$Iso_type),
        mae      = as.numeric(.data$mae),
        model    = model_name
      )
  }

  model_label <- c(
    overall  = "Overall",
    serospec = "Sero-specific",
    combined = "Combined flexneri"
  )

  ov  <- standardize_mae(mae_overall,  "overall")
  sp  <- standardize_mae(mae_serospec, "serospec")
  cmb <- standardize_mae(mae_combined, "combined")

  summarize_model <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::group_by(.data$antigen, .data$Iso_type) |>
      dplyr::summarise(
        med = stats::median(.data$mae, na.rm = TRUE),
        q25 = stats::quantile(.data$mae, 0.25, na.rm = TRUE),
        q75 = stats::quantile(.data$mae, 0.75, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        value = sprintf("%.2f (%.2f\u2013%.2f)", .data$med, .data$q25, .data$q75)
      ) |>
      dplyr::select("antigen", "Iso_type", "value")
  }

  ov_sum  <- summarize_model(ov)
  sp_sum  <- summarize_model(sp)
  cmb_sum <- summarize_model(cmb)

  base_rows <- dplyr::bind_rows(
    if (!is.null(ov_sum))  ov_sum  |> dplyr::select("antigen", "Iso_type"),
    if (!is.null(sp_sum))  sp_sum  |> dplyr::select("antigen", "Iso_type"),
    if (!is.null(cmb_sum)) cmb_sum |> dplyr::select("antigen", "Iso_type")
  ) |>
    dplyr::distinct()

  safe_join <- function(tbl, x, col_name) {
    if (is.null(x)) {
      return(dplyr::mutate(tbl, !!col_name := NA_character_))
    }
    dplyr::left_join(tbl, dplyr::rename(x, !!col_name := "value"),
                     by = c("antigen", "Iso_type"))
  }

  tbl <- base_rows |>
    safe_join(ov_sum,  "Overall") |>
    safe_join(sp_sum,  "Sero-specific") |>
    safe_join(cmb_sum, "Combined flexneri")

  if (is.null(n_lookup)) {
    n_lookup <- dplyr::bind_rows(ov, sp, cmb) |>
      dplyr::group_by(.data$antigen) |>
      dplyr::summarise(n = dplyr::n_distinct(.data$sid), .groups = "drop")
  } else {
    n_lookup <- tibble::as_tibble(n_lookup) |>
      dplyr::transmute(antigen = as.character(.data$antigen),
                       n       = as.integer(.data$n))
  }
  tbl <- dplyr::left_join(tbl, n_lookup, by = "antigen")

  mae_long <- dplyr::bind_rows(ov, sp, cmb)

  calc_best_model <- function(dat) {
    available_models <- sort(unique(dat$model))
    if (length(available_models) == 1) {
      n_ids <- dplyr::n_distinct(dat$sid)
      return(tibble::tibble(
        best_display = sprintf("%s (%d)", model_label[available_models[1]], n_ids)
      ))
    }
    ids_by_model <- dat |>
      dplyr::group_by(.data$model) |>
      dplyr::summarise(ids = list(sort(unique(.data$sid))), .groups = "drop")
    common_ids <- Reduce(intersect, ids_by_model$ids)
    if (length(common_ids) == 0) {
      return(tibble::tibble(best_display = "N/A"))
    }
    comp_dat <- dat |> dplyr::filter(.data$sid %in% common_ids)
    winners <- comp_dat |>
      dplyr::group_by(.data$sid) |>
      dplyr::slice_min(.data$mae, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    win_count <- winners |>
      dplyr::count(.data$model, name = "wins") |>
      dplyr::arrange(dplyr::desc(.data$wins))
    best_key  <- win_count$model[1]
    best_wins <- win_count$wins[1]
    tibble::tibble(
      best_display = sprintf("%s (%d/%d)",
                             model_label[best_key], best_wins,
                             length(common_ids))
    )
  }

  best_tbl <- mae_long |>
    dplyr::group_by(.data$antigen, .data$Iso_type) |>
    dplyr::group_modify(~ calc_best_model(.x)) |>
    dplyr::ungroup() |>
    dplyr::rename(`Best model` = "best_display")

  tbl <- dplyr::left_join(tbl, best_tbl, by = c("antigen", "Iso_type")) |>
    dplyr::mutate(
      dplyr::across(
        c("Overall", "Sero-specific", "Combined flexneri", "Best model"),
        ~ ifelse(is.na(.), "N/A", .)
      )
    )

  display_df <- tbl |>
    dplyr::mutate(
      antigen  = factor(.data$antigen,
                        levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")),
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA"))
    ) |>
    dplyr::arrange(.data$antigen, .data$Iso_type) |>
    dplyr::select(
      Antigen  = "antigen",
      Isotype  = "Iso_type",
      "n",
      "Overall",
      "Sero-specific",
      "Combined flexneri",
      "Best model"
    )

  flextable::flextable(display_df) |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 7, part = "all") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1:2, align = "left", part = "body") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption(
      caption = paste(
        "Per-individual prediction accuracy (MAE) across modeling approaches.",
        "Values are median MAE (IQR) on the log-MFI scale.",
        "Best model determined using shared individuals across compared models."
      )
    )
}


#' Build a Flextable of Kinetic Parameters with Formatted Headers
#'
#' Converts a display-ready parameter tibble (from [format_param_table()])
#' to a `flextable` with properly subscripted/superscripted column headers
#' and Greek-letter symbols.
#'
#' @param display_df A tibble as returned by [format_param_table()], with
#'   columns `Biomarker`, `y0`, `y1`, `t1`, `alpha`, `rho`.
#' @param caption Character scalar table caption.
#' @param footer_lines Character vector of footer annotation lines.
#'
#' @return A `flextable` object.
#'
#' @examples
#' \dontrun{
#' tbl <- prep_pop_params(overall_IpaB_pop_6, "IpaB") |>
#'   transform_pop_params() |>
#'   summarise_pop_params() |>
#'   format_param_table() |>
#'   build_kinetic_flextable(
#'     caption = "IpaB kinetic parameters.",
#'     footer_lines = "Values are posterior medians (95% CrI)."
#'   )
#' }
#'
#' @importFrom flextable flextable theme_booktabs fontsize bold align
#'   set_table_properties set_caption add_footer_lines padding
#'   line_spacing compose as_paragraph as_chunk as_sub as_sup
#' @export
build_kinetic_flextable <- function(display_df,
                                    caption      = "",
                                    footer_lines = character()) {

  ft <- flextable::flextable(display_df) |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 7, part = "body") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "body") |>
    flextable::align(align = "center", part = "header") |>
    flextable::align(j = 1, align = "left", part = "body") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption(caption = caption)

  if (length(footer_lines) > 0) {
    ft <- flextable::add_footer_lines(ft, values = footer_lines) |>
      flextable::fontsize(size = 6, part = "footer") |>
      flextable::line_spacing(part = "footer", space = 0.7) |>
      flextable::padding(
        part = "footer", padding.top = 0, padding.bottom = 0
      )
  }

  # subscripted/superscripted column headers
  if ("y0" %in% names(display_df)) {
    ft <- flextable::compose(
      ft, part = "header", j = "y0",
      value = flextable::as_paragraph(
        flextable::as_chunk("y"), flextable::as_sub("0")
      )
    )
  }
  if ("y1" %in% names(display_df)) {
    ft <- flextable::compose(
      ft, part = "header", j = "y1",
      value = flextable::as_paragraph(
        flextable::as_chunk("y"), flextable::as_sub("1")
      )
    )
  }
  if ("t1" %in% names(display_df)) {
    ft <- flextable::compose(
      ft, part = "header", j = "t1",
      value = flextable::as_paragraph(
        flextable::as_chunk("t"), flextable::as_sub("1")
      )
    )
  }
  if ("alpha" %in% names(display_df)) {
    ft <- flextable::compose(
      ft, part = "header", j = "alpha",
      value = flextable::as_paragraph(
        flextable::as_chunk("\u03B1"),
        flextable::as_chunk(" (year"),
        flextable::as_sup("-1"),
        flextable::as_chunk(")")
      )
    )
  }
  if ("rho" %in% names(display_df)) {
    ft <- flextable::compose(
      ft, part = "header", j = "rho",
      value = flextable::as_paragraph(
        flextable::as_chunk("\u03C1")
      )
    )
  }

  ft
}
