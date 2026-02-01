#' Internal utility functions
#'
#' @keywords internal
#' @name utils_internal
NULL

#' Get time variable name from a case_data object
#'
#' @param dataset A serodynamics case_data object
#' @return Character scalar with the time variable name
#' @keywords internal
#' @noRd
get_timeindays_var <- function(dataset) {
  time_attrs <- c("timeindays", "time_in_days", "time")
  time_var <- NULL
  
  # Check attributes first
  for (attr_name in time_attrs) {
    attr_val <- attr(dataset, attr_name, exact = TRUE)
    if (!is.null(attr_val)) {
      time_var <- attr_val
      break
    }
  }
  
  # If not found in attributes, look for common column names
  if (is.null(time_var)) {
    col_names <- names(dataset)
    for (possible_name in time_attrs) {
      if (possible_name %in% col_names) {
        time_var <- possible_name
        break
      }
    }
  }
  
  if (is.null(time_var)) {
    stop("Could not determine time variable from dataset")
  }
  
  time_var
}

#' Antibody kinetics model function
#'
#' Evaluates the antibody trajectory model at specified time points.
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
  # Antibody kinetics model:
  # - Linear rise from y0 to y1 over [0, t1]
  # - Weibull-like decay after t1
  
  result <- numeric(length(t))
  
  # Rising phase (t <= t1)
  rising <- t <= t1
  result[rising] <- y0 + (y1 - y0) * (t[rising] / t1)
  
  # Decay phase (t > t1)
  decaying <- t > t1
  result[decaying] <- y0 + (y1 - y0) * exp(-alpha * ((t[decaying] - t1)^shape))
  
  result
}
