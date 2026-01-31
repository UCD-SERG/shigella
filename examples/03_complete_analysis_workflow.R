# Example 3: Complete Analysis Workflow
# This script demonstrates a complete end-to-end analysis workflow
# combining data processing, incidence estimation, and visualization

library(shigella)
library(dplyr)
library(tibble)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load and Process Data from Multiple Regions
# -----------------------------------------------------------------------------

set.seed(789)

# Create mock multi-region dataset
mock_data <- tibble(
  study_name = rep("SOSAR", 300),
  site_name = sample(c("MA USA", "Ghana", "Niger"), 300, replace = TRUE),
  sid = rep(1:100, 3),
  isotype_name = sample(c("IgG", "IgA"), 300, replace = TRUE),
  age = runif(300, 0.5, 60),
  n_ipab_MFI = rlnorm(300, log(1000), 1)
)

# Extract cross-sectional data for each region
regions <- c("MA USA", "Ghana", "Niger")
xs_data_list <- list()

for (region in regions) {
  xs_data_list[[region]] <- create_xs_data(
    df = mock_data,
    filter_countries = region,
    filter_antigen_iso = c("IgG"),
    value_col = n_ipab_MFI
  )
}

cat("Data extracted for", length(regions), "regions\n")

# -----------------------------------------------------------------------------
# 2. Summary Statistics by Region
# -----------------------------------------------------------------------------

summary_by_region <- lapply(names(xs_data_list), function(region) {
  xs_data_list[[region]] %>%
    summarise(
      Region = region,
      n = n(),
      mean_age = mean(age),
      mean_value = mean(value),
      median_value = median(value),
      sd_value = sd(value)
    )
}) %>%
  bind_rows()

cat("\nSummary statistics by region:\n")
print(summary_by_region)

# -----------------------------------------------------------------------------
# 3. Age Distribution by Region
# -----------------------------------------------------------------------------

age_dist_plot <- bind_rows(xs_data_list, .id = "Region") %>%
  ggplot(aes(x = age, fill = Region)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity") +
  facet_wrap(~Region, ncol = 1) +
  labs(
    title = "Age Distribution by Region",
    x = "Age (years)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

cat("\nAge distribution plot created\n")
print(age_dist_plot)

# -----------------------------------------------------------------------------
# 4. Antibody Level Distribution by Age Category
# -----------------------------------------------------------------------------

antibody_plot <- bind_rows(xs_data_list, .id = "Region") %>%
  ggplot(aes(x = ageCat, y = value, fill = Region)) +
  geom_boxplot() +
  facet_wrap(~Region, ncol = 3) +
  scale_y_log10() +
  labs(
    title = "Antibody Levels by Age Category and Region",
    x = "Age Category",
    y = "Antibody Level (MFI, log scale)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cat("Antibody level distribution plot created\n")
print(antibody_plot)

# -----------------------------------------------------------------------------
# 5. Mock Incidence Estimation (Placeholder)
# -----------------------------------------------------------------------------
# In a real analysis, you would:
# 1. Fit antibody decay curves to longitudinal data
# 2. Use those curves to estimate incidence from cross-sectional data
# For now, we create mock estimates

mock_estimates <- list(
  USA = list(incidence.rate = 0.12, SE = 0.03, CI.lwr = 0.06, CI.upr = 0.18),
  Ghana = list(incidence.rate = 0.25, SE = 0.05, CI.lwr = 0.15, CI.upr = 0.35),
  Niger = list(incidence.rate = 0.18, SE = 0.04, CI.lwr = 0.10, CI.upr = 0.26)
)

# Create incidence table (would normally use create_incidence_table() with real estimates)
incidence_table <- tibble(
  Country = names(mock_estimates),
  Incidence_Rate = sapply(mock_estimates, function(x) x$incidence.rate),
  SE = sapply(mock_estimates, function(x) x$SE),
  CI.lwr = sapply(mock_estimates, function(x) x$CI.lwr),
  CI.upr = sapply(mock_estimates, function(x) x$CI.upr)
)

cat("\nEstimated incidence rates by region:\n")
print(incidence_table)

# -----------------------------------------------------------------------------
# 6. Visualize Incidence Estimates
# -----------------------------------------------------------------------------

incidence_plot <- ggplot(incidence_table, aes(x = Country, y = Incidence_Rate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(
    aes(ymin = CI.lwr, ymax = CI.upr),
    width = 0.2,
    color = "steelblue"
  ) +
  labs(
    title = "Estimated Shigella Incidence Rates by Region",
    x = "Region",
    y = "Incidence Rate (per person-year)",
    caption = "Error bars represent 95% confidence intervals"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0, face = "italic")
  )

cat("\nIncidence rate plot created\n")
print(incidence_plot)

# -----------------------------------------------------------------------------
# 7. Summary Metrics Table for Sample Size Analysis
# -----------------------------------------------------------------------------

# Create mock summary metrics (would come from simulation studies)
mock_metrics <- expand.grid(
  sample_size = c(100, 200, 300, 400),
  Age_Group = c("<5", "5-15", "16+")
) %>%
  as_tibble() %>%
  mutate(
    empirical_se = 0.1 / sqrt(sample_size) * runif(n(), 0.8, 1.2)
  )

# Create plot
se_plot <- plot_summary_metrics(mock_metrics)

cat("\nSample size vs. empirical SE plot created\n")
print(se_plot)

# -----------------------------------------------------------------------------
# 8. Export Results (Optional)
# -----------------------------------------------------------------------------

# In a real analysis, you might save results:
# write.csv(incidence_table, "incidence_estimates.csv", row.names = FALSE)
# write.csv(summary_by_region, "regional_summaries.csv", row.names = FALSE)
# ggsave("incidence_plot.png", incidence_plot, width = 8, height = 6)

cat("\nAnalysis workflow complete!\n")
cat("\nNext steps with real data:\n")
cat("1. Load real Shigella data from data/ directory\n")
cat("2. Fit antibody decay curves using JAGS/Stan\n")
cat("3. Estimate incidence using serocalculator/serodynamics\n")
cat("4. Run power analysis simulations\n")
cat("5. Compare estimates across regions and age groups\n")
