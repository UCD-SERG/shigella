#' Internal utility functions
#'
#' @keywords internal
#' @name utils_internal
NULL

# NOTE: ab() and get_timeindays_var() are internal re-implementations of
# functions from the serodynamics package (ucdavis/serodynamics).
# They are duplicated here to avoid using ::: (non-exported access),
# which is not allowed in R packages passing R CMD check.
# If serodynamics updates these functions, this file must be updated to match.

#' Get time variable name from a case_data object
#'
#' @param dataset A serodynamics case_data object
#' @return Character scalar with the time variable name
#' @keywords internal
#' @noRd
get_timeindays_var <- function(dataset) {
  var <- attr(dataset, "timeindays")
  if (is.null(var)) var <- "timeindays"
  return(var)
}

#' Compute log-linear rise rate for antibody kinetics
#'
#' Helper that computes the exponential rise rate beta = log(y1/y0) / t1,
#' used internally by `ab()`.
#'
#' @param y0 Baseline antibody level
#' @param y1 Peak antibody level
#' @param t1 Time to peak (days)
#' @return Numeric scalar: the log-linear rise rate
#' @keywords internal
#' @noRd
bt <- function(y0, y1, t1) {
  to_return <- log(y1 / y0) / t1
  return(to_return)
}

#' Antibody kinetics trajectory function
#'
#' Evaluates the antibody trajectory model at specified time points.
#' Mirrors `serodynamics:::ab()` -- see ucdavis/serodynamics/R/ab.R.
#'
#' @param t Numeric vector of time points (days)
#' @param y0 Baseline antibody level
#' @param y1 Peak antibody level
#' @param t1 Time to peak (days)
#' @param alpha Decay rate parameter
#' @param shape Decay shape parameter (rho)
#' @return Numeric vector of predicted antibody levels
#' @keywords internal
#' @noRd
ab <- function(t, y0, y1, t1, alpha, shape) {
  beta <- bt(y0, y1, t1)
  yt <- ifelse(
    t <= t1,
    y0 * exp(beta * t),
    (y1^(1 - shape) - (1 - shape) * alpha * (t - t1))^(1 / (1 - shape))
  )
  return(yt)
}
