#' Extract "newperson" parameter draws and summarize for Table 2
#'
#' @param draws_long Posterior draws in long format with columns:
#'   Subject, Iso_type, Chain, Iteration, Parameter, value.
#' @param antigen_label Character scalar antigen name to attach.
#'
#' @return A tibble of newperson draws in wide format with columns:
#' Iteration, Chain, antigen, Iso_type, y0, y1, t1, alpha, rho.
#'
#' @export
prep_newperson_params <- function(draws_long, antigen_label) {
  draws_long |>
    dplyr::filter(.data$Subject == "newperson") |>
    dplyr::mutate(antigen = antigen_label) |>
    dplyr::select(.data$Iteration, .data$Chain, .data$antigen,
                  .data$Iso_type, .data$Parameter, .data$value) |>
    tidyr::pivot_wider(names_from = .data$Parameter, values_from = .data$value) |>
    # The model output uses "shape" for the decay-shape parameter, which is renamed to "rho" for consistency
    dplyr::rename(rho = .data$shape) |>
    dplyr::select(.data$Iteration, .data$Chain, .data$antigen, .data$Iso_type,
                  .data$y0, .data$y1, .data$t1, .data$alpha, .data$rho)
}
