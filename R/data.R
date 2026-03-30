# Mock data documentation

#' Mock posterior draws for testing
#'
#' A mock dataset of posterior parameter draws in long format, mimicking the
#' structure expected by functions like \code{\link{compute_residual_metrics}}
#' and \code{\link{predict_posterior_at_times}}.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Subject}{Character. Subject ID (e.g., "newperson", "SOSAR-22008")}
#'   \item{Iso_type}{Character. Isotype (e.g., "IgG", "IgA")}
#'   \item{Chain}{Integer. MCMC chain number}
#'   \item{Iteration}{Integer. MCMC iteration number}
#'   \item{Parameter}{Character. Parameter name (y0, y1, t1, alpha, shape)}
#'   \item{value}{Numeric. Parameter value}
#' }
#'
#' @details
#' This is synthetic data generated for testing and examples. Real Shigella
#' posterior draws will be added to the package separately.
#'
#' Parameters represent:
#' \itemize{
#'   \item \code{y0}: Baseline antibody level
#'   \item \code{y1}: Peak antibody level
#'   \item \code{t1}: Time to peak (days)
#'   \item \code{alpha}: Decay rate parameter
#'   \item \code{shape}: Decay shape parameter (rho)
#' }
#'
#' @examples
#' head(mock_posterior_draws)
#' table(mock_posterior_draws$Parameter)
"mock_posterior_draws"

#' Mock case data for testing
#'
#' A mock longitudinal antibody dataset compatible with
#' \code{serodynamics::as_case_data()}, for testing functions like
#' \code{\link{compute_residual_metrics}}.
#'
#' @format A data frame with class \code{c("case_data", "data.frame")}
#'   and columns:
#' \describe{
#'   \item{id}{Character. Subject ID}
#'   \item{antigen_iso}{Character. Isotype (e.g., "IgG", "IgA")}
#'   \item{timepoint}{Numeric. Time in days since infection}
#'   \item{value}{Numeric. Antibody measurement}
#' }
#'
#' @details
#' This is synthetic data generated for testing and examples. Real Shigella
#' case data will be added separately.
#'
#' The dataset has attributes:
#' \itemize{
#'   \item \code{attr(mock_case_data, "timeindays") = "timepoint"}
#'   \item \code{attr(mock_case_data, "value_var") = "value"}
#' }
#'
#' @examples
#' head(mock_case_data)
#' table(mock_case_data$id)
"mock_case_data"
