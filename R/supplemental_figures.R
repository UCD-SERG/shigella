#' Age-Stratified Spaghetti Plot of Raw Antibody Responses (Supplemental)
#'
#' Produces a faceted spaghetti plot of raw longitudinal IgG and IgA antibody
#' trajectories for one age subgroup (either children <5 or ≥5 years),
#' coloured by infecting serotype.  Used for S1 Fig (under-5) and S2 Fig
#' (plus-5) in the Chapter 1 manuscript.
#'
#' @param dL_list A named list of `case_data` objects, one per antigen,
#'   in the order `IpaB`, `Sf2a`, `Sf3a`, `Sf6`, `Sonnei`.  Each object
#'   must have `sid` (or `id`), `timeindays`, `result`, `isotype_name`
#'   (or `antigen_iso`), and `cohort_name` columns.
#' @param age_label Character scalar used as the plot subtitle, e.g.,
#'   `"<5 years"` or `"\u22655 years"`.
#' @param col_map Named character vector of hex colours for infecting
#'   serotypes.  Default matches the main spaghetti figure.
#' @param alpha_map Named numeric vector of line alpha values per serotype.
#'   Default matches the main spaghetti figure.
#' @param lw_map Named numeric vector of line widths per serotype.
#'   Default matches the main spaghetti figure.
#'
#' @return A `ggplot` object with `facet_grid(antigen ~ isotype_name)`.
#'
#' @details
#' The function calls `.bind_antigen_list()` internally to combine the named
#' list into a single long tibble before plotting.
#'
#' @seealso [fig5_with_individuals()], [build_figure4_row()]
#'
#' @examples
#' \dontrun{
#' data("dL_clean_Ipab_new_under5",   package = "shigella")
#' data("dL_clean_sf2a_new_under5",   package = "shigella")
#' data("dL_clean_sf3a_new_under5",   package = "shigella")
#' data("dL_clean_sf6_new_under5",    package = "shigella")
#' data("dL_clean_sonnei_new_under5", package = "shigella")
#'
#' dL_under5 <- list(
#'   IpaB   = dL_clean_Ipab_new_under5,
#'   Sf2a   = dL_clean_sf2a_new_under5,
#'   Sf3a   = dL_clean_sf3a_new_under5,
#'   Sf6    = dL_clean_sf6_new_under5,
#'   Sonnei = dL_clean_sonnei_new_under5
#' )
#' plot_age_spaghetti(dL_under5, age_label = "<5 years")
#' }
#'
#' @importFrom dplyr bind_rows mutate case_when filter
#' @importFrom purrr imap_dfr
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual
#'   scale_linewidth_manual scale_alpha_manual scale_y_log10 facet_grid
#'   labs theme_bw theme element_rect element_text element_blank
#' @export
plot_age_spaghetti <- function(dL_list,
                                age_label = "",
                                col_map = c(
                                  Sf2a   = "#D81B60",
                                  sonnei = "#1E88E5",
                                  Sf3a   = "#E65100",
                                  Sf6    = "#2E7D32",
                                  Other  = "grey45"
                                ),
                                alpha_map = c(
                                  Sf2a   = 0.45,
                                  sonnei = 0.40,
                                  Sf3a   = 0.50,
                                  Sf6    = 0.45,
                                  Other  = 0.20
                                ),
                                lw_map = c(
                                  Sf2a   = 0.5,
                                  sonnei = 0.4,
                                  Sf3a   = 0.45,
                                  Sf6    = 0.40,
                                  Other  = 0.2
                                )) {

  antigen_levels <- c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")

  dL_all <- purrr::imap_dfr(dL_list, function(dat, antigen_name) {
    .extract_individual_trajectories(dat) |>
      dplyr::mutate(
        antigen      = antigen_name,
        cohort_name  = if ("cohort_name" %in% names(dat))
          dat$cohort_name[match(.data$id, dat[[.id_col(dat)]])]
        else NA_character_
      )
  }) |>
    dplyr::mutate(
      antigen = factor(.data$antigen, levels = antigen_levels),
      isotype_name = factor(.data$iso, levels = c("IgG", "IgA")),
      infecting_serotype = dplyr::case_when(
        .data$cohort_name == "Sf2a"   ~ "Sf2a",
        .data$cohort_name == "sonnei" ~ "sonnei",
        .data$cohort_name == "Sf3a"   ~ "Sf3a",
        .data$cohort_name == "Sf6"    ~ "Sf6",
        TRUE                          ~ "Other"
      ) |>
        factor(levels = c("Sf2a", "sonnei", "Sf3a", "Sf6", "Other"))
    ) |>
    dplyr::filter(!is.na(.data$t), !is.na(.data$value))

  ggplot2::ggplot(
    dL_all,
    ggplot2::aes(
      x     = .data$t,
      y     = .data$value,
      group = .data$id,
      color = .data$infecting_serotype
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        linewidth = .data$infecting_serotype,
        alpha     = .data$infecting_serotype
      ),
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = .data$infecting_serotype),
      size = 1.2, shape = 16
    ) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::scale_linewidth_manual(values = lw_map, guide = "none") +
    ggplot2::scale_alpha_manual(values = alpha_map, guide = "none") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(antigen ~ isotype_name, scales = "free_y") +
    ggplot2::labs(
      x        = "Days since symptom onset",
      y        = "Antibody level (Normalized MFI, log scale)",
      color    = "Infecting serotype",
      subtitle = age_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text       = ggplot2::element_text(colour = "white", face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom",
      plot.subtitle    = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}


#' Per-Individual MAE Slopegraph Across Modelling Approaches (S4 Fig)
#'
#' Produces a faceted slopegraph where each line connects one individual's
#' MAE across the overall, serotype-specific, and (where applicable) combined
#' *S. flexneri* models.  Green lines indicate that the alternative model
#' reduced prediction error; grey lines indicate the overall model performed
#' better.  Corresponds to S4 Fig in the Chapter 1 manuscript.
#'
#' @param mae_overall A tibble of per-individual MAE from the overall model,
#'   as returned by [get_mae()].
#' @param mae_serospec A tibble of per-individual MAE from serotype-specific
#'   models, or `NULL`.
#' @param mae_combined A tibble of per-individual MAE from the combined
#'   *S. flexneri* model, or `NULL`.
#'
#' @return A `ggplot` object with `facet_wrap(~ biomarker)`, two panels per
#'   row.  Only participants shared across all compared models within each
#'   panel are plotted.
#'
#' @examples
#' \dontrun{
#' p_slope <- plot_mae_slopegraph(
#'   mae_overall  = mae_overall,
#'   mae_serospec = mae_serospec,
#'   mae_combined = mae_combined
#' )
#' print(p_slope)
#' }
#'
#' @importFrom dplyr bind_rows mutate filter group_by group_modify ungroup
#'   summarise
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual
#'   scale_fill_manual scale_shape_manual facet_wrap labs theme_minimal
#'   theme element_text element_blank element_rect
#' @export
plot_mae_slopegraph <- function(mae_overall,
                                 mae_serospec = NULL,
                                 mae_combined = NULL) {

  mae_long <- dplyr::bind_rows(
    mae_overall  |> dplyr::mutate(model = "Overall"),
    if (!is.null(mae_serospec))
      mae_serospec |> dplyr::mutate(model = "Sero-specific"),
    if (!is.null(mae_combined))
      mae_combined |> dplyr::mutate(model = "Combined")
  ) |>
    dplyr::mutate(
      biomarker = paste0(.data$antigen, " ", .data$Iso_type),
      model     = factor(.data$model,
                         levels = c("Overall", "Sero-specific", "Combined"))
    )

  # Keep only antigen-isotypes with >1 model
  multi_biomarkers <- mae_long |>
    dplyr::group_by(.data$biomarker) |>
    dplyr::summarise(
      n_models = dplyr::n_distinct(.data$model),
      .groups  = "drop"
    ) |>
    dplyr::filter(.data$n_models > 1) |>
    dplyr::pull(.data$biomarker)

  mae_slope_filtered <- mae_long |>
    dplyr::filter(.data$biomarker %in% multi_biomarkers)

  # Restrict to individuals shared across all models within each biomarker
  mae_shared <- mae_slope_filtered |>
    dplyr::group_by(.data$biomarker) |>
    dplyr::group_modify(~ {
      ids_by_model <- .x |>
        dplyr::group_by(.data$model) |>
        dplyr::summarise(ids = list(unique(.data$sid)), .groups = "drop")
      common_ids <- Reduce(intersect, ids_by_model$ids)
      .x |> dplyr::filter(.data$sid %in% common_ids)
    }) |>
    dplyr::ungroup()

  # Determine whether the alternative beats the overall per individual
  mae_direction <- mae_shared |>
    dplyr::group_by(.data$biomarker, .data$sid) |>
    dplyr::mutate(
      overall_mae  = .data$mae[.data$model == "Overall"][1],
      best_alt_mae = min(.data$mae[.data$model != "Overall"], na.rm = TRUE),
      improved     = dplyr::if_else(
        .data$best_alt_mae < .data$overall_mae,
        "Alternative better",
        "Overall better"
      )
    ) |>
    dplyr::ungroup()

  ggplot2::ggplot(
    mae_direction,
    ggplot2::aes(x = .data$model, y = .data$mae, group = .data$sid)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(color = .data$improved),
      alpha     = 0.6,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$model, fill = .data$model),
      size   = 2.4,
      stroke = 0.2,
      color  = "black",
      alpha  = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Overall better"      = "grey70",
        "Alternative better"  = "#5BBE7A"
      ),
      name = "Lower MAE"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Overall"       = "#3B82F6",
        "Sero-specific" = "#EF4444",
        "Combined"      = "#10B981"
      ),
      name = "Model"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Overall"       = 22L,
        "Sero-specific" = 21L,
        "Combined"      = 24L
      ),
      name = "Model"
    ) +
    ggplot2::facet_wrap(~ biomarker, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title    = "Per-individual MAE across modelling approaches",
      subtitle = paste(
        "Each line connects the same individual across available models.",
        "Only individuals shared across all compared models are shown."
      ),
      x = "Modelling approach",
      y = "Per-individual MAE (log-MFI scale)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(
        face = "bold", size = 13, hjust = 0.5
      ),
      plot.subtitle    = ggplot2::element_text(size = 10, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text       = ggplot2::element_text(
        colour = "white", face = "bold"
      ),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.box         = "vertical",
      axis.text.x        = ggplot2::element_text(angle = 20, hjust = 1)
    )
}


# ── Internal helpers ─────────────────────────────────────────────────────────

#' @keywords internal
.id_col <- function(dat) {
  if ("id"  %in% names(dat)) "id"  else
    if ("sid" %in% names(dat)) "sid" else "id"
}
