#' Posterior antibody predictions at a set of times for given subjects
#'
#' Evaluates the two-phase curve serodynamics:::ab() for every posterior draw of
#' the named subjects at the requested times.
#'
#' @param model A fitted `sr_model` (individual-level draws).
#' @param ids Character vector of subject ids.
#' @param antigen_iso Isotype to evaluate (`"IgG"`/`"IgA"`).
#' @param times Numeric vector of times (days).
#' @return Long tibble of per-draw predictions in column `res`.
#' @export
predict_posterior_at_times <- function(model, ids, antigen_iso, times) {
  sr_model_sub <- model |>
    dplyr::filter(.data$Subject %in% ids, .data$Iso_type == antigen_iso)

  param_medians_wide <- sr_model_sub |>
    dplyr::select(dplyr::all_of(c("Chain", "Iteration", "Iso_type",
                                  "Parameter", "value", "Subject"))) |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
    dplyr::arrange(.data$Chain, .data$Iteration) |>
    dplyr::mutate(
      antigen_iso = factor(.data$Iso_type),
      id = as.factor(.data$Subject),
      r = .data$shape
    ) |>
    dplyr::select(-dplyr::all_of(c("Iso_type", "Subject")))

  if (!"sample_id" %in% names(param_medians_wide)) {
    param_medians_wide <- dplyr::mutate(param_medians_wide,
                                        sample_id = dplyr::row_number())
  }

  dt1 <- data.frame(t = times) |>
    dplyr::mutate(idx = dplyr::row_number()) |>
    tidyr::pivot_wider(names_from = "idx", values_from = "t",
                       names_prefix = "time") |>
    dplyr::slice(rep(seq_len(dplyr::n()), each = nrow(param_medians_wide)))

  cbind(param_medians_wide, dt1) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("time"), values_to = "t") |>
    dplyr::select(-dplyr::all_of("name")) |>
    # TODO: replace with exported serodynamics API once ab() is exported (tracked separately).
    dplyr::mutate(res = serodynamics:::ab(.data$t, .data$y0, .data$y1,
                                          .data$t1, .data$alpha, .data$shape))
}
