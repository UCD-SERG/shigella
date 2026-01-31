#' Create Noise Parameter Data Frame for Serodynamics Models
#'
#' Generates a tibble containing biological and measurement noise parameters
#' for a specific geographic region/country, used in seroincidence estimation.
#'
#' @param country Character string specifying the country or geographic region
#'   name (e.g., "MA USA", "Ghana", "Niger").
#'
#' @return A tibble with one row containing noise parameters:
#'   \item{antigen_iso}{Antigen-isotype combination (default: "IgG")}
#'   \item{Country}{Factor with the specified country name}
#'   \item{y.low}{Lower measurement bound (default: 25)}
#'   \item{eps}{Measurement error rate (default: 0.25)}
#'   \item{nu}{Biological noise parameter (default: 0.5)}
#'   \item{y.high}{Upper measurement bound (default: 200000)}
#'
#' @examples
#' \dontrun{
#' # Create noise parameters for USA
#' noise_usa <- create_noise_df("MA USA")
#'
#' # Create noise parameters for Ghana
#' noise_ghana <- create_noise_df("Ghana")
#' }
#'
#' @importFrom tibble tibble
#' @export
create_noise_df <- function(country) {
  noise_df <- tibble(
    antigen_iso = c("IgG"),
    Country     = factor(c(country)), # Use input argument for Country
    y.low       = c(25),
    eps         = c(0.25),
    nu          = c(0.5),
    y.high      = c(200000)
  )
  
  return(noise_df)
}