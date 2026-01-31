#' Plot Empirical Standard Error vs Sample Size
#'
#' Creates a line plot showing the relationship between sample size and
#' empirical standard error, with separate lines for different age groups.
#' Useful for visualizing simulation results and power analysis.
#'
#' @param summary_metrics A data frame containing summary statistics from
#'   simulation studies, with columns:
#'   \itemize{
#'     \item `sample_size` - Sample size values
#'     \item `empirical_se` - Empirical standard error estimates
#'     \item `Age_Group` - Age group categories for stratification
#'   }
#'
#' @return A ggplot object displaying empirical standard error vs. sample size,
#'   colored by age group.
#'
#' @examples
#' \dontrun{
#' # Create mock summary metrics
#' mock_metrics <- data.frame(
#'   sample_size = rep(c(100, 200, 300, 400), each = 3),
#'   empirical_se = runif(12, 0.01, 0.1),
#'   Age_Group = rep(c("<5", "5-15", "16+"), 4)
#' )
#'
#' # Generate plot
#' plot_summary_metrics(mock_metrics)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @export
plot_summary_metrics <- function(summary_metrics) {
  summary_metrics |> 
    ggplot() +
    aes(x = sample_size, y = empirical_se, color = Age_Group) +
    geom_line() +
    geom_point() +
    labs(
      title = "Empirical Standard Error vs. Sample Size",
      x = "Sample Size",
      y = "Empirical Standard Error",
      color = "Age Group"
    ) +
    theme_minimal()
}