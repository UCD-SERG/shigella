# Internal helpers

# Coerce one get_mae() tibble to a uniform sid/antigen/Iso_type/mae/model frame
# (NULL passes through; missing required columns abort with a clear message).
#' @keywords internal
#' @noRd
.mae_standardize <- function(df, model_name) {
  required_cols <- c("sid", "antigen", "Iso_type", "mae")
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- tibble::as_tibble(df)
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) cli::cli_abort(paste0(
    "Input '", model_name, "' missing: ", paste(miss, collapse = ", ")))
  dplyr::transmute(df, sid = as.character(.data$sid),
                   antigen = as.character(.data$antigen),
                   Iso_type = as.character(.data$Iso_type),
                   mae = as.numeric(.data$mae), model = model_name)
}

# Median (IQR) MAE string per antigen x isotype for one model (NULL -> NULL).
#' @keywords internal
#' @noRd
.mae_summarise_one <- function(df) {
  if (is.null(df)) return(NULL)
  df |>
    dplyr::group_by(.data$antigen, .data$Iso_type) |>
    dplyr::summarise(value = sprintf("%.2f (%.2f\u2013%.2f)",
      stats::median(.data$mae, na.rm = TRUE),
      stats::quantile(.data$mae, 0.25, na.rm = TRUE),
      stats::quantile(.data$mae, 0.75, na.rm = TRUE)), .groups = "drop")
}

# Assemble the wide median(IQR) table: each present model becomes one column;
# absent models get an all-NA column so the display step stays uniform.
#' @keywords internal
#' @noRd
.mae_wide <- function(sums) {
  present <- sums[!vapply(sums, is.null, logical(1))]
  base <- purrr::map(present, ~ dplyr::select(.x, "antigen", "Iso_type")) |>
    dplyr::bind_rows() |>
    dplyr::distinct()
  for (nm in names(sums)) {
    if (is.null(sums[[nm]])) {
      base[[nm]] <- NA_character_
    } else {
      base <- dplyr::left_join(base, dplyr::rename(sums[[nm]], !!nm := "value"),
                               by = c("antigen", "Iso_type"))
    }
  }
  base
}

# For one antigen x isotype group, the model that wins on the most shared
# individuals (with win count and number shared).
#' @keywords internal
#' @noRd
.mae_best_one <- function(dat, model_label) {
  models <- sort(unique(dat$model))
  if (length(models) == 1)
    return(tibble::tibble(best_display = model_label[models[1]],
                          wins = NA_integer_, n_shared = NA_integer_))
  ids_by <- dat |> dplyr::group_by(.data$model) |>
    dplyr::summarise(ids = list(sort(unique(.data$sid))), .groups = "drop")
  common <- Reduce(intersect, ids_by$ids)
  if (length(common) == 0)
    return(tibble::tibble(best_display = "N/A", wins = NA_integer_, n_shared = NA_integer_)) # nolint: line_length_linter.
  winners <- dat |> dplyr::filter(.data$sid %in% common) |>
    dplyr::group_by(.data$sid) |>
    dplyr::slice_min(.data$mae, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::count(.data$model, name = "wins") |>
    dplyr::arrange(dplyr::desc(.data$wins))
  tibble::tibble(best_display = model_label[winners$model[1]],
                 wins = winners$wins[1], n_shared = length(common))
}

# Per-row "antigen-iso: best for k/N" string (Results text / S5 Table).
#' @keywords internal
#' @noRd
.mae_winner_footnote <- function(best_tbl) {
  best_tbl |>
    dplyr::filter(!is.na(.data$wins)) |>
    dplyr::transmute(label = sprintf("%s-%s: %s best for %d/%d",
      .data$antigen, .data$Iso_type, .data$best_display, .data$wins, .data$n_shared)) |> # nolint: line_length_linter.
    dplyr::pull(.data$label) |>
    paste(collapse = "; ")
}

# Join the best-model column, fill N/A, order rows, select display columns.
#' @keywords internal
#' @noRd
.mae_display_df <- function(wide_tbl, best_tbl) {
  wide_tbl |>
    dplyr::left_join(dplyr::select(best_tbl, "antigen", "Iso_type",
                                   `Best model` = "best_display"),
                     by = c("antigen", "Iso_type")) |>
    dplyr::mutate(dplyr::across(c("Overall", "Serotype-specific",
                                  "Combined flexneri", "Best model"),
                                ~ dplyr::if_else(is.na(.), "N/A", .))) |>
    dplyr::mutate(antigen = factor(.data$antigen,
                                   levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")), # nolint: line_length_linter.
                  Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA"))) |> # nolint: line_length_linter.
    dplyr::arrange(.data$antigen, .data$Iso_type) |>
    dplyr::select(Antigen = "antigen", Isotype = "Iso_type",
                  "Overall", "Serotype-specific", "Combined flexneri", "Best model") # nolint: line_length_linter.
}

# flextable styling for the MAE table.
#' @keywords internal
#' @noRd
.mae_flextable <- function(display_df) {
  flextable::flextable(display_df) |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 7, part = "all") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1:2, align = "left", part = "body") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption(paste0(
      "Per-individual prediction accuracy (MAE) across modeling approaches. ",
      "Values are median MAE (IQR) on the log-MFI scale. The 'Best model' column ", # nolint: line_length_linter.
      "identifies the modeling approach with the lowest median MAE among ",
      "individuals shared across compared models."))
}

# Main function

#' Table 4: per-individual prediction accuracy (MAE) across modeling approaches
#'
#' Orchestrates the Table 4 helpers: standardise each model's MAE
#' (`.mae_standardize()`), summarise to median (IQR) (`.mae_summarise_one()`),
#' assemble the wide table (`.mae_wide()`), pick the best model per
#' antigen-isotype (`.mae_best_one()`), then build the display frame
#' (`.mae_display_df()`) and flextable (`.mae_flextable()`).
#'
#' @param mae_overall,mae_serospec,mae_combined `get_mae()` tibbles per model
#'   class (`NULL` if not applicable).
#' @return A list with `flextable` (the rendered table) and `winner_footnote`
#'   (a per-row "best for k/N" string for the Results text / S5 Table).
#' @export
table4_mae <- function(mae_overall, mae_serospec = NULL, mae_combined = NULL) {
  model_label <- c(overall = "Overall", serospec = "Serotype-specific",
                   combined = "Combined flexneri")

  ov  <- .mae_standardize(mae_overall,  "overall")
  sp  <- .mae_standardize(mae_serospec, "serospec")
  cmb <- .mae_standardize(mae_combined, "combined")

  wide <- .mae_wide(list(
    Overall             = .mae_summarise_one(ov),
    `Serotype-specific` = .mae_summarise_one(sp),
    `Combined flexneri` = .mae_summarise_one(cmb)
  ))

  best_tbl <- dplyr::bind_rows(ov, sp, cmb) |>
    dplyr::group_by(.data$antigen, .data$Iso_type) |>
    dplyr::group_modify(~ .mae_best_one(.x, model_label)) |>
    dplyr::ungroup()

  display_df <- .mae_display_df(wide, best_tbl)

  list(flextable       = .mae_flextable(display_df),
       winner_footnote = .mae_winner_footnote(best_tbl))
}
