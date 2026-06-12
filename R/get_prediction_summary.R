#' Posterior predictive summary (median + 95% CrI) over a time grid
#'
#' Used by the individual-trajectory figures to draw a fitted curve with a
#' credible band for one subject x isotype.
#'
#' @param model_output Fitted `sr_model` (individual-level draws).
#' @param sid Subject id.
#' @param iso Isotype (`"IgG"`/`"IgA"`).
#' @param times Time grid (days).
#' @return Tibble `t, med, lo, hi`.
#' @export
get_prediction_summary <- function(model_output, sid, iso,
                                    times = seq(0, 200, by = 1)) {
  sub <- model_output |>
    dplyr::filter(.data$Subject == sid, .data$Iso_type == iso)

  if (nrow(sub) == 0) {
    cli::cli_warn("No data for {.val {sid}} {iso}")
    return(tibble::tibble(t = numeric(), med = numeric(),
                          lo = numeric(), hi = numeric()))
  }

  wide <- sub |>
    dplyr::select("Iteration", "Chain", "Parameter", "value") |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value")

  # TODO: replace with exported serodynamics API once ab() is exported
  # (tracked separately).
  pred_mat <- vapply(
    times,
    function(tt) serodynamics:::ab(tt, wide$y0, wide$y1, wide$t1,
                                   wide$alpha, wide$shape),
    numeric(nrow(wide))
  )

  tibble::tibble(
    t   = times,
    med = apply(pred_mat, 2, stats::median, na.rm = TRUE),
    lo  = apply(pred_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    hi  = apply(pred_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}
