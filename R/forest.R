#' Forest Plot of Population-Level Kinetic Parameters Across Modelling Approaches
#'
#' Produces a multi-panel forest plot (one panel per kinetic parameter)
#' comparing posterior medians and 95% credible intervals from the overall,
#' serotype-specific, and (optionally) combined *S. flexneri* hierarchical
#' models.  Corresponds to S6 Fig in the Chapter 1 manuscript.
#'
#' @param pop_overall A summarised parameter tibble (from
#'   [summarise_pop_params()]) for the overall model (n = 48), containing all
#'   five antigens.
#' @param pop_serospec A summarised parameter tibble for the serotype-specific
#'   models (*S. flexneri* 2a, 3a, and *S. sonnei*).
#' @param pop_combined A summarised parameter tibble for the combined
#'   *S. flexneri* model (n = 25), or `NULL` (default) to omit this group.
#' @param colors A named character vector of hex colors for the three model
#'   groups.  Names must match `"Overall"`, `"Serotype-specific"`, and
#'   `"Combined flexneri"`.
#' @param point_size Numeric.  Size of the posterior median point.
#'   Default `3`.
#' @param line_width Numeric.  Width of the credible-interval line.
#'   Default `0.8`.
#' @param dodge_width Numeric.  Vertical dodge between model groups on each
#'   antigen row.  Default `0.5`.
#'
#' @return A `patchwork` object combining five parameter panels stacked
#'   vertically, with a shared legend at the bottom of the last panel.
#'
#' @details
#' Each panel corresponds to one kinetic parameter:
#' \describe{
#'   \item{Panel 1}{Baseline antibody level \eqn{y_0}}
#'   \item{Panel 2}{Peak antibody level \eqn{y_1}}
#'   \item{Panel 3}{Time to peak \eqn{t_1} (days)}
#'   \item{Panel 4}{Decay rate \eqn{\alpha} (year\eqn{^{-1}})}
#'   \item{Panel 5}{Decay shape \eqn{\rho}}
#' }
#' Parameters \eqn{y_0}, \eqn{y_1}, \eqn{t_1}, and \eqn{\alpha} are
#' plotted on a log10 scale; \eqn{\rho} is plotted on a linear scale.
#' Within each panel, rows represent antigens and points are colour-coded
#' by model group.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(patchwork)
#'
#' pop_ov <- bind_rows(
#'   prep_pop_params(overall_IpaB_pop_6,   "IpaB"),
#'   prep_pop_params(overall_Sf2a_pop_6,   "Sf2a"),
#'   prep_pop_params(overall_Sf3a_pop_6,   "Sf3a"),
#'   prep_pop_params(overall_Sf6_pop_6,    "Sf6"),
#'   prep_pop_params(overall_Sonnei_pop_6, "Sonnei")
#' ) |> transform_pop_params() |> summarise_pop_params()
#'
#' pop_sp <- bind_rows(
#'   prep_pop_params(serotype_sf2a_3,   "Sf2a"),
#'   prep_pop_params(serotype_sf3a_3,   "Sf3a"),
#'   prep_pop_params(serotype_sonnei_3, "Sonnei")
#' ) |> transform_pop_params() |> summarise_pop_params()
#'
#' pop_cmb <- bind_rows(
#'   prep_pop_params(combined_flexneri_sf2a, "Sf2a"),
#'   prep_pop_params(combined_flexneri_sf3a, "Sf3a")
#' ) |> transform_pop_params() |> summarise_pop_params()
#'
#' p_forest <- plot_model_comparison_forest(
#'   pop_overall  = pop_ov,
#'   pop_serospec = pop_sp,
#'   pop_combined = pop_cmb
#' )
#' print(p_forest)
#' }
#'
#' @importFrom dplyr bind_rows mutate select group_by group_split
#' @importFrom purrr map
#' @importFrom ggplot2 ggplot aes geom_linerange geom_point
#'   scale_colour_manual scale_x_log10 facet_wrap labs theme_minimal
#'   theme element_blank element_text element_rect margin
#'   position_dodge
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
plot_model_comparison_forest <- function(pop_overall,
                                          pop_serospec,
                                          pop_combined = NULL,
                                          colors = c(
                                            "Overall"           = "#2563EB",
                                            "Serotype-specific" = "#DC2626",
                                            "Combined flexneri" = "#16A34A"
                                          ),
                                          point_size   = 3,
                                          line_width   = 0.8,
                                          dodge_width  = 0.5) {

  # --- Reshape each model's summary to long format ------------------------
  reshape_to_long <- function(pop_table, model_label) {
    dplyr::bind_rows(
      pop_table |>
        dplyr::select("antigen", "Iso_type",
                      med = "y0_med", lo = "y0_lo", hi = "y0_hi") |>
        dplyr::mutate(param = "y[0]~(baseline)", param_order = 1L),
      pop_table |>
        dplyr::select("antigen", "Iso_type",
                      med = "y1_med", lo = "y1_lo", hi = "y1_hi") |>
        dplyr::mutate(param = "y[1]~(peak)", param_order = 2L),
      pop_table |>
        dplyr::select("antigen", "Iso_type",
                      med = "t1_med", lo = "t1_lo", hi = "t1_hi") |>
        dplyr::mutate(param = "t[1]~(days~to~peak)", param_order = 3L),
      pop_table |>
        dplyr::select("antigen", "Iso_type",
                      med = "alpha_med", lo = "alpha_lo", hi = "alpha_hi") |>
        dplyr::mutate(param = "alpha~(decay~rate)", param_order = 4L),
      pop_table |>
        dplyr::select("antigen", "Iso_type",
                      med = "rho_med", lo = "rho_lo", hi = "rho_hi") |>
        dplyr::mutate(param = "rho~(shape)", param_order = 5L)
    ) |>
      dplyr::mutate(model = model_label)
  }

  params_long <- dplyr::bind_rows(
    reshape_to_long(pop_overall,  "Overall"),
    reshape_to_long(pop_serospec, "Serotype-specific"),
    if (!is.null(pop_combined))
      reshape_to_long(pop_combined, "Combined flexneri")
    else
      NULL
  ) |>
    dplyr::mutate(
      antigen  = factor(.data$antigen,
                        levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")),
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA")),
      model    = factor(.data$model,
                        levels = c("Overall",
                                   "Serotype-specific",
                                   "Combined flexneri")),
      param    = reorder(.data$param, .data$param_order)
    )

  log_params <- 1:4   # y0, y1, t1, alpha on log scale

  make_panel <- function(df_sub, use_log = TRUE) {
    p <- ggplot2::ggplot(
      df_sub,
      ggplot2::aes(x = .data$med, y = .data$antigen, colour = .data$model)
    ) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
        linewidth = line_width,
        position  = ggplot2::position_dodge(width = dodge_width),
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        size     = point_size,
        position = ggplot2::position_dodge(width = dodge_width)
      ) +
      ggplot2::scale_colour_manual(values = colors, name = NULL) +
      ggplot2::facet_wrap(~ Iso_type, ncol = 2) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        strip.background   = ggplot2::element_rect(fill = "grey20"),
        strip.text         = ggplot2::element_text(
          colour = "white", face = "bold", size = 12
        ),
        axis.text.y        = ggplot2::element_text(face = "bold", size = 10),
        legend.position    = "none",
        plot.margin        = ggplot2::margin(4, 8, 4, 4)
      )

    if (use_log) {
      p <- p + ggplot2::scale_x_log10(
        labels = function(x) {
          ifelse(x >= 1,
                 formatC(x, format = "g", digits = 3),
                 formatC(x, format = "e", digits = 0))
        }
      )
    }
    p
  }

  panels <- params_long |>
    dplyr::group_by(.data$param, .data$param_order) |>
    dplyr::group_split() |>
    purrr::map(function(df_sub) {
      ord     <- unique(df_sub$param_order)
      use_log <- ord %in% log_params
      make_panel(df_sub, use_log = use_log) +
        ggplot2::labs(
          title = parse(text = unique(as.character(df_sub$param)))
        )
    })

  # Legend on the last panel only
  panels[[length(panels)]] <- panels[[length(panels)]] +
    ggplot2::theme(legend.position = "bottom")

  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title = paste(
        "Population-level antibody kinetic parameters:",
        "Overall vs. Serotype-specific models"
      ),
      subtitle = paste(
        "Posterior median and 95% credible interval.",
        "Blue = overall (n=48); Red = serotype-specific;",
        "Green = combined S. flexneri (n=25, where applicable)"
      ),
      caption = paste(
        "Panels faceted by isotype (IgG | IgA).\n",
        "Parameters y\u2080, y\u2081, t\u2081, and \u03B1 on log scale.\n",
        "IpaB and Sf6: overall model only.\n",
        "MCMC: 4 chains, 15,000 retained samples per chain."
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(
          face = "bold", size = 15, hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 9, hjust = 0.5, colour = "grey40"
        ),
        plot.caption  = ggplot2::element_text(
          size = 10, colour = "grey50", hjust = 0.5
        )
      )
    )
}
