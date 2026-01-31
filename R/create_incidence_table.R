#' Create Incidence Rate Summary Table from Multiple Estimates
#'
#' Extracts and combines estimated incidence rates from multiple seroincidence
#' estimates (typically from different geographic regions) into a single tidy
#' summary table.
#'
#' @param ... Named seroincidence estimate objects (output from
#'   `serocalculator::est.incidence()` or `serodynamics::est.incidence()`).
#'   Names will be used as country/region identifiers in the output table.
#'
#' @return A tibble with two columns:
#'   \item{Country}{Character vector of country/region names (from input names)}
#'   \item{Incidence_Rate}{Numeric vector of estimated incidence rates}
#'
#' @examples
#' \dontrun{
#' # Assuming estimates exist for multiple regions
#' incidence_summary <- create_incidence_table(
#'   USA = est_usa,
#'   Ghana = est_ghana,
#'   Niger = est_niger
#' )
#' }
#'
#' @importFrom tibble tibble
#' @export
# create table of incidence.rate of each region
create_incidence_table <- function(...) {
  # Capture input objects and their names
  est_list <- list(...)
  country_names <- names(est_list)
  
  # Extract incidence.rate from the summary() of each estimate object
  incidence_rates <- sapply(est_list, function(x) summary(x)$incidence.rate)
  
  # Create a tidy tibble
  incidence_table <- tibble(
    Country = country_names,
    Incidence_Rate = incidence_rates
  )
  
  return(incidence_table)
}