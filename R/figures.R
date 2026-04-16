#' Build One Row of the Individual-Level Model Comparison Figure
#'
#' Produces a two-panel ggplot (IgG | IgA) showing predicted antibody
#' trajectories under two or three modeling approaches for a single
#' representative participant, with observed data overlaid.
#'
#' @param sid A single character subject ID.
#' @param model_overall Fitted overall model object (n = 48).
#' @param model_sero Fitted serotype-specific model object.
#' @param data_overall Case_data object for the overall model pool.
#' @param data_sero Case_data object for the serotype-specific pool.
#' @param model_combined Fitted combined *S. flexneri* model object, or
#'   `NULL` (default) for a two-way comparison.
#' @param row_label A character string used as the row title.
#' @param antigen_label Short antigen name displayed in panel titles.
#' @param xlim Numeric vector of length 2 for x-axis limits. Default
#'   `c(0, 200)`.
#' @param times Numeric vector of prediction times. Default
#'   `seq(0, 200, by = 1)`.
#'
#' @return A `patchwork` object combining IgG and IgA panels.
#'
#' @examples
#' \dontrun{
#' row_sf2a <- build_figure4_row(
#'   sid            = "SOSAR-22008",
#'   model_overall  = overall_Sf2a_pop_6,
#'   model_sero     = serotype_sf2a_3,
#'   data_overall   = dL_clean_sf2a,
#'   data_sero      = dL_serotype_sf2a,
#'   row_label      = "A) S. flexneri 2a (n = 17)",
#'   antigen_label  = "SF2a"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   scale_color_manual scale_linetype_manual scale_y_log10
#'   coord_cartesian labs theme_minimal theme element_text element_blank
#'   element_rect
#' @importFrom scales label_comma
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
build_figure4_row <- function(sid,
                               model_overall,
                               model_sero,
                               data_overall,
                               data_sero,
                               model_combined = NULL,
                               row_label      = "",
                               antigen_label  = "Sf2a",
                               xlim           = c(0, 200),
                               times          = seq(0, 200, by = 1)) {

  col_overall  <- "#2563EB"
  col_sero     <- "#DC2626"
  col_combined <- "#16A34A"

  make_panel <- function(iso) {
    obs <- get_observed(data_sero, sid, iso)
    if (nrow(obs) == 0) obs <- get_observed(data_overall, sid, iso)

    pred_ovr  <- get_prediction_summary(model_overall, sid, iso, times)
    pred_sero <- get_prediction_summary(model_sero,    sid, iso, times)

    panel_title <- paste0(antigen_label, " \u2013 ", iso, " \u2013 ", sid)

    color_values <- c(
      "Overall"           = col_overall,
      "Sero-specific"     = col_sero,
      "Combined flexneri" = col_combined
    )
    lt_values <- c(
      "Overall"           = "solid",
      "Sero-specific"     = "dashed",
      "Combined flexneri" = "dotted"
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = pred_ovr,
        ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
        fill = col_overall, alpha = 0.12
      ) +
      ggplot2::geom_line(
        data = pred_ovr,
        ggplot2::aes(x = .data$t, y = .data$med,
                     color = "Overall", linetype = "Overall"),
        linewidth = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = pred_sero,
        ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
        fill = col_sero, alpha = 0.12
      ) +
      ggplot2::geom_line(
        data = pred_sero,
        ggplot2::aes(x = .data$t, y = .data$med,
                     color = "Sero-specific", linetype = "Sero-specific"),
        linewidth = 0.9
      )

    if (!is.null(model_combined)) {
      pred_cmb <- get_prediction_summary(model_combined, sid, iso, times)
      if (nrow(pred_cmb) > 0) {
        p <- p +
          ggplot2::geom_ribbon(
            data = pred_cmb,
            ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
            fill = col_combined, alpha = 0.10
          ) +
          ggplot2::geom_line(
            data = pred_cmb,
            ggplot2::aes(x = .data$t, y = .data$med,
                         color = "Combined flexneri",
                         linetype = "Combined flexneri"),
            linewidth = 0.9
          )
      }
    }

    if (nrow(obs) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = obs,
          ggplot2::aes(x = .data$t, y = .data$value),
          size = 2.8, shape = 21, fill = "white",
          color = "black", stroke = 0.8
        )
    }

    p +
      ggplot2::scale_color_manual(
        name = "Model", values = color_values,
        breaks = names(color_values)
      ) +
      ggplot2::scale_linetype_manual(
        name = "Model", values = lt_values,
        breaks = names(lt_values)
      ) +
      ggplot2::scale_y_log10(labels = scales::label_comma()) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(
        x = "Days since symptom onset",
        y = "MFI (log scale)",
        title = panel_title
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(
          face = "bold", size = 10, hjust = 0.5
        ),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey20"),
        strip.text       = ggplot2::element_text(
          colour = "white", face = "bold"
        ),
        legend.position  = "none"
      )
  }

  p_igg <- make_panel("IgG")
  p_iga <- make_panel("IgA")

  (p_igg | p_iga) +
    patchwork::plot_annotation(
      title = row_label,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0)
      )
    )
}


#' Plot Population-Level Kinetic Trajectories with Individual Observations
#'
#' Creates a two-facet (IgG | IgA) ggplot showing the population-level
#' posterior median trajectory and 95% credible interval, overlaid with
#' individual observed antibody trajectories in gray.
#'
#' @param model_obj A fitted model object with a `"population_params"`
#'   attribute.
#' @param raw_dataset A case_data object for individual trajectories.
#' @param antigen_label Character scalar for the panel title.
#' @param t_grid Numeric vector of times. Default `seq(0, 210, by = 5)`.
#' @param line_color Hex color for population line. Default `"#1f77b4"`.
#' @param ribbon_alpha Alpha for CrI ribbon. Default `0.15`.
#' @param individual_alpha Alpha for individual lines. Default `0.12`.
#' @param individual_color Color for individual lines. Default `"grey70"`.
#' @param log_y Log10 y-axis? Default `TRUE`.
#' @param xlim X-axis limits. Default `c(0, 210)`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' fig5_with_individuals(
#'   serotype_sonnei_3, dL_serotype_sonnei,
#'   "B) S. sonnei (Sero-specific, n=11)"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon facet_wrap labs
#'   theme_minimal theme element_text scale_y_log10 coord_cartesian
#' @importFrom scales label_comma
#' @importFrom dplyr filter mutate summarise
#' @importFrom tidyr pivot_wider crossing
#' @importFrom purrr map_dfr
#' @export
fig5_with_individuals <- function(model_obj,
                                   raw_dataset,
                                   antigen_label,
                                   t_grid            = seq(0, 210, by = 5),
                                   line_color        = "#1f77b4",
                                   ribbon_alpha      = 0.15,
                                   individual_alpha  = 0.12,
                                   individual_color  = "grey70",
                                   log_y             = TRUE,
                                   xlim              = c(0, 210)) {

  pop_curves <- .extract_pop_curves(model_obj, t_grid)
  indiv_data <- .extract_individual_trajectories(raw_dataset)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = indiv_data,
      ggplot2::aes(x = .data$t, y = .data$value, group = .data$id),
      color = individual_color, alpha = individual_alpha, linewidth = 0.3
    ) +
    ggplot2::geom_ribbon(
      data = pop_curves,
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
      fill = "grey50", alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(
      data = pop_curves,
      ggplot2::aes(x = .data$t, y = .data$med),
      color = line_color, linewidth = 1.0, lineend = "round"
    ) +
    ggplot2::facet_wrap(~ iso, ncol = 2) +
    ggplot2::labs(
      x = "Days since symptom onset",
      y = "Normalized MFI",
      title = antigen_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text  = ggplot2::element_text(face = "bold", size = 12),
      axis.title  = ggplot2::element_text(size = 11),
      axis.text   = ggplot2::element_text(size = 10),
      legend.position = "none"
    )

  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_comma())
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}


#' Plot IpaB Kinetics with Age-Stratified Population Curves
#'
#' Overlays three population trajectories: overall (blue solid), <5 years
#' (orange dashed), and >=5 years (green dashed).
#'
#' @param model_overall Fitted overall IpaB model (n = 48).
#' @param model_under5 Fitted IpaB model for children <5 years.
#' @param model_plus5 Fitted IpaB model for children >=5 years.
#' @param raw_overall Case_data for all participants.
#' @param raw_under5 Case_data for children <5 years.
#' @param raw_plus5 Case_data for children >=5 years.
#' @param antigen_label Panel title. Default
#'   `"A) IpaB (Overall + age-stratified)"`.
#' @param t_grid Prediction times. Default `seq(0, 210, by = 5)`.
#' @param log_y Log10 y-axis? Default `TRUE`.
#' @param xlim X-axis limits. Default `c(0, 210)`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' fig5_ipab_with_age(
#'   overall_IpaB_pop_6,
#'   overall_IpaB_pop_under5,
#'   overall_IpaB_pop_plus5,
#'   dL_clean_Ipab,
#'   dL_clean_Ipab_new_under5,
#'   dL_clean_Ipab_new_plus5
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon facet_wrap labs
#'   theme_minimal theme element_text scale_color_manual
#'   scale_linetype_manual scale_linewidth_manual scale_y_log10
#'   coord_cartesian guides guide_legend unit
#' @importFrom scales label_comma
#' @importFrom dplyr filter mutate bind_rows
#' @export
fig5_ipab_with_age <- function(model_overall,
                                model_under5,
                                model_plus5,
                                raw_overall,
                                raw_under5,
                                raw_plus5,
                                antigen_label = "A) IpaB (Overall + age-stratified)",
                                t_grid        = seq(0, 210, by = 5),
                                log_y         = TRUE,
                                xlim          = c(0, 210)) {

  curve_overall <- .extract_pop_curves(model_overall, t_grid) |>
    dplyr::mutate(group = "Overall (n=48)")
  curve_under5  <- .extract_pop_curves(model_under5, t_grid) |>
    dplyr::mutate(group = "<5 years")
  curve_plus5   <- .extract_pop_curves(model_plus5, t_grid) |>
    dplyr::mutate(group = "\u22655 years")

  all_curves <- dplyr::bind_rows(curve_overall, curve_under5, curve_plus5) |>
    dplyr::mutate(
      group = factor(.data$group,
                     levels = c("<5 years", "\u22655 years", "Overall (n=48)"))
    )

  indiv_data <- .extract_individual_trajectories(raw_overall)

  color_vals <- c(
    "<5 years"       = "#E64A19",
    "\u22655 years"  = "#2E7D32",
    "Overall (n=48)" = "#1f77b4"
  )
  lt_vals <- c(
    "<5 years"       = "dashed",
    "\u22655 years"  = "dashed",
    "Overall (n=48)" = "solid"
  )
  alpha_ribbon <- c(
    "<5 years"       = 0.10,
    "\u22655 years"  = 0.10,
    "Overall (n=48)" = 0.08
  )

  p <- ggplot2::ggplot()

  p <- p +
    ggplot2::geom_line(
      data = indiv_data,
      ggplot2::aes(x = .data$t, y = .data$value, group = .data$id),
      color = "grey75", alpha = 0.18, linewidth = 0.3
    )

  for (grp in c("Overall (n=48)", "<5 years", "\u22655 years")) {
    grp_data <- all_curves |> dplyr::filter(.data$group == grp)
    p <- p +
      ggplot2::geom_ribbon(
        data = grp_data,
        ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
        fill  = color_vals[grp],
        alpha = alpha_ribbon[grp]
      )
  }

  p <- p +
    ggplot2::geom_line(
      data = all_curves,
      ggplot2::aes(
        x = .data$t, y = .data$med,
        color = .data$group, linewidth = .data$group, linetype = .data$group
      ),
      lineend = "round"
    ) +
    ggplot2::scale_color_manual(values = color_vals, name = "Age group") +
    ggplot2::scale_linewidth_manual(
      values = c("<5 years" = 1.0, "\u22655 years" = 1.0,
                 "Overall (n=48)" = 1.0),
      guide = "none"
    ) +
    ggplot2::scale_linetype_manual(values = lt_vals, name = "Age group") +
    ggplot2::guides(
      color    = ggplot2::guide_legend(
        override.aes = list(linewidth = 1.2)
      ),
      linetype = ggplot2::guide_legend()
    ) +
    ggplot2::facet_wrap(~ iso, ncol = 2) +
    ggplot2::labs(
      x = "Days since symptom onset", y = "Normalized MFI",
      title = antigen_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(
        face = "bold", size = 13, hjust = 0.5
      ),
      strip.text       = ggplot2::element_text(face = "bold", size = 12),
      axis.title       = ggplot2::element_text(size = 11),
      axis.text        = ggplot2::element_text(size = 10),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = 10),
      legend.text      = ggplot2::element_text(size = 9),
      legend.key.width = ggplot2::unit(1.5, "cm")
    )

  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_comma())
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}


# ── Internal helpers ─────────────────────────────────────────────────────────

#' Extract Population-Level Curves from a Model Object (Internal)
#'
#' Computes posterior median and 95% CrI for the population-level
#' antibody trajectory using the `.ab()` internal function.
#'
#' @param model_obj Fitted model with `"population_params"` attribute.
#' @param t_grid Numeric vector of time points.
#' @return A tibble with columns `Iso_type`, `t`, `med`, `lo`, `hi`, `iso`.
#' @keywords internal
#' @noRd
.extract_pop_curves <- function(model_obj, t_grid) {
  pop <- attr(model_obj, "population_params")
  isotypes <- c("IgG", "IgA")

  purrr::map_dfr(isotypes, function(iso) {
    mu_wide <- pop |>
      dplyr::filter(
        .data$Population_Parameter == "mu.par",
        .data$Iso_type == iso
      ) |>
      dplyr::select("Chain", "Iteration", "Iso_type", "Parameter", "value") |>
      tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
      dplyr::mutate(
        y0_nat    = exp(.data$y0),
        y1_nat    = exp(.data$y0) + exp(.data$y1),
        t1_nat    = exp(.data$t1),
        alpha_nat = exp(.data$alpha),
        shape_nat = exp(.data$shape) + 1
      )

    mu_wide |>
      tidyr::crossing(t = t_grid) |>
      dplyr::mutate(
        res = .ab(
          .data$t, .data$y0, .data$y1,
          .data$t1, .data$alpha, .data$shape
        )
      ) |>
      dplyr::filter(is.finite(.data$res), .data$res > 0) |>
      dplyr::summarise(
        med = stats::quantile(.data$res, 0.50,  na.rm = TRUE),
        lo  = stats::quantile(.data$res, 0.025, na.rm = TRUE),
        hi  = stats::quantile(.data$res, 0.975, na.rm = TRUE),
        .by = c("Iso_type", "t")
      )
  }) |>
    dplyr::mutate(iso = factor(.data$Iso_type, levels = c("IgG", "IgA")))
}


#' Extract Individual Trajectories from a Case-Data Object (Internal)
#' @keywords internal
#' @noRd
.extract_individual_trajectories <- function(raw_dataset) {
  time_col  <- if ("timeindays"    %in% names(raw_dataset)) "timeindays"    else
    if ("Actual day"    %in% names(raw_dataset)) "Actual day"    else "timepoint"
  value_col <- if ("result"        %in% names(raw_dataset)) "result"        else {
    mfi <- names(raw_dataset)[grepl("_MFI$", names(raw_dataset))]
    if (length(mfi) == 1) mfi else "result"
  }
  id_col    <- if ("sid"           %in% names(raw_dataset)) "sid"           else "id"
  iso_col   <- if ("isotype_name"  %in% names(raw_dataset)) "isotype_name"  else
    if ("Iso_type"      %in% names(raw_dataset)) "Iso_type"      else "antigen_iso"

  raw_dataset |>
    dplyr::transmute(
      id    = .data[[id_col]],
      t     = .data[[time_col]],
      value = .data[[value_col]],
      iso   = .data[[iso_col]]
    ) |>
    dplyr::filter(!is.na(.data$value), !is.na(.data$t)) |>
    dplyr::mutate(iso = factor(.data$iso, levels = c("IgG", "IgA")))
}
