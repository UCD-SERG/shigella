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