#' Posterior predictions at specified times for given subjects and antigen/isotype
#'
#' Generates draw-level posterior predictions of the antibody trajectory at user-specified
#' time points, for one or more subjects and a selected antigen/isotype. This is a low-level
#' helper used by residual-based posterior predictive diagnostics.
#'
#' @param model A data frame of posterior draws in long format with columns:
#'   \code{Subject}, \code{Iso_type}, \code{Chain}, \code{Iteration}, \code{Parameter}, \code{value}.
#' @param ids Character vector of subject IDs to include (matched against \code{Subject}).
#' @param antigen_iso Character scalar specifying the antigen/isotype to include
#'   (matched against \code{Iso_type}).
#' @param times Numeric vector of time points (days) at which to evaluate predictions.
#'
#' @return A tibble with one row per (posterior draw \eqn{\times} time \eqn{\times} subject),
#' including the evaluated prediction \code{res}. Output includes at least:
#' \describe{
#'   \item{id}{Subject ID (character).}
#'   \item{t}{Time (days) at which prediction was evaluated.}
#'   \item{Chain}{MCMC chain index (if present in \code{model}).}
#'   \item{Iteration}{MCMC iteration index (if present in \code{model}).}
#'   \item{sample_id}{Row index for the draw (added if missing).}
#'   \item{y0, y1, t1, alpha, shape}{Model parameters (wide).}
#'   \item{res}{Predicted antibody level at time \code{t}.}
#' }
#'
#' @details
#' This function pivots posterior draws to wide format (parameters as columns),
#' expands them over \code{times}, and evaluates the antibody curve via
#' an internal implementation of the antibody kinetics model using parameters
#' \code{y0}, \code{y1}, \code{t1}, \code{alpha}, and \code{shape}.
#'
#' @seealso \code{\link{compute_residual_metrics}}
#'
#' @examples
#' \dontrun{
#' preds <- predict_posterior_at_times(
#'   model = overall_sf2a,
#'   ids = "newperson",
#'   antigen_iso = "IgG",
#'   times = c(0, 30, 90, 180)
#' )
#' }
#'
#' @keywords internal
predict_posterior_at_times <- function(model, ids, antigen_iso, times) {

  sr_model_sub <- model |>
    dplyr::filter(.data$Subject %in% ids, .data$Iso_type == antigen_iso)

  param_wide <- sr_model_sub |>
    dplyr::select(.data$Chain, .data$Iteration, .data$Iso_type,
                  .data$Parameter, .data$value, .data$Subject) |>
    tidyr::pivot_wider(names_from = .data$Parameter, values_from = .data$value) |>
    dplyr::arrange(.data$Chain, .data$Iteration) |>
    dplyr::mutate(
      antigen_iso = factor(.data$Iso_type),
      id = as.character(.data$Subject)
    ) |>
    dplyr::select(-c("Iso_type", "Subject"))

  if (!"sample_id" %in% names(param_wide)) {
    param_wide <- param_wide |>
      dplyr::mutate(sample_id = dplyr::row_number())
  }

  dt <- tibble::tibble(t = times) |>
    dplyr::mutate(idx = dplyr::row_number()) |>
    tidyr::pivot_wider(names_from = "idx", values_from = "t", names_prefix = "time") |>
    dplyr::slice(rep(seq_len(dplyr::n()), each = nrow(param_wide)))

  predictions <- cbind(param_wide, dt) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("time"), values_to = "t") |>
    dplyr::select(-"name") |>
    dplyr::mutate(
      res = ab(.data$t, .data$y0, .data$y1, .data$t1, .data$alpha, .data$shape)
    )

  predictions
}
