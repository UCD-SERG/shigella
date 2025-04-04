#' Calculate Empirical Standard Error for Incidence Rates
#'
#' This function computes the empirical standard error (SE) of incidence rates
#' from the provided dataset and adds sample size and age group information.
#'
#' @param data A data frame containing a column named \code{incidence.rate}.
#' @param sample_size Integer representing the sample size used in the simulation.
#' @param age_group Character string specifying the age group (e.g., "Age 0-2").
#'
#' @return A data frame with \code{sample_size}, \code{empirical_se}, and \code{Age_Group}.
#' @examples
#' \dontrun{
#' data <- data.frame(incidence.rate = c(0.1, 0.2, 0.15, 0.18))
#' calculate_metrics(data, sample_size = 100, age_group = "Age 0-2")
#' }
#' @export
calculate_metrics <- function(data, sample_size = attr(data, "sample_size, age_group = attr(data, "age_group") {
  data.frame(
    sample_size = sample_size,
    empirical_se = sd(data$incidence.rate, na.rm = TRUE),
    Age_Group = age_group
  )
}
