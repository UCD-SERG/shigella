# ==========================================================================
# 03_analyze_results.R
#
# Analyze Phase 2 simulation results.
# Produces:
#   - Summary table: bias, RMSE, coverage per scenario
#   - Recovery plots: true vs estimated rho_B
#   - Divergent transitions report
# ==========================================================================

setwd("~/chapter2")

library(dplyr)
library(tidyr)
library(ggplot2)

# ==========================================================================
# 1. Load simulation results
# ==========================================================================
all_results <- readRDS("outputs/02_simulation_results.rds")

# Flatten to data frame
results_df <- dplyr::bind_rows(
  lapply(all_results, function(scn_list) {
    dplyr::bind_rows(scn_list)
  })
) |>
  dplyr::filter(status == "OK")

cat("Total successful fits:", nrow(results_df), "\n")
cat("Failed fits:", sum(unlist(lapply(all_results, function(s)
  sum(sapply(s, function(r) r$status == "FAILED"))))), "\n\n")

# Check available columns
cat("Columns in results_df:\n")
print(names(results_df))

# Add n_divergent if it was not saved in the simulation results
if (!"n_divergent" %in% names(results_df)) {
  warning("Column 'n_divergent' not found in results_df. Setting n_divergent = 0.")
  results_df <- results_df |>
    dplyr::mutate(n_divergent = 0)
}

# ==========================================================================
# 2. Summary metrics per scenario
# ==========================================================================
summary_metrics <- results_df |>
  dplyr::group_by(scenario, true_rho_B) |>
  dplyr::summarise(
    n_reps = dplyr::n(),
    mean_estimate = mean(est_rho_B_median),
    bias = mean(bias),
    rmse = sqrt(mean(bias^2)),
    coverage_95 = mean((true_rho_B >= est_rho_B_lo) &
                       (true_rho_B <= est_rho_B_hi)),
    mean_ci_width = mean(est_rho_B_hi - est_rho_B_lo),
    total_divergent = sum(n_divergent),
    pct_divergent = mean(n_divergent > 0) * 100,
    median_runtime_min = median(elapsed_min),
    .groups = "drop"
  )

cat("=== Summary metrics ===\n")
print(summary_metrics)

saveRDS(summary_metrics, "outputs/03_summary_metrics.rds")
write.csv(summary_metrics, "outputs/03_summary_metrics.csv", row.names = FALSE)

# ==========================================================================
# 3. Recovery plot
# ==========================================================================
p_recovery <- ggplot(results_df,
                     aes(x = true_rho_B, y = est_rho_B_median,
                         color = scenario)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey50") +
  geom_errorbar(aes(ymin = est_rho_B_lo, ymax = est_rho_B_hi),
                width = 0.02, alpha = 0.4) +
  geom_jitter(width = 0.01, height = 0, size = 2.5, alpha = 0.7) +
  scale_color_manual(values = c("A" = "#2166AC",
                                 "B" = "#92C5DE",
                                 "C" = "#B2182B")) +
  labs(
    title = "Phase 2: Sigma_B Parameter Recovery",
    subtitle = "Each point = 1 simulation replicate, 95% CrI as error bars",
    x = "True rho_B",
    y = "Estimated rho_B (posterior median)"
  ) +
  theme_bw(base_size = 12) +
  facet_wrap(~ scenario, scales = "fixed")

ggsave("outputs/03_recovery_plot.png", p_recovery,
       width = 12, height = 5, dpi = 150, bg = "white")

# ==========================================================================
# 4. Divergent transitions report
# ==========================================================================
div_report <- results_df |>
  dplyr::group_by(scenario) |>
  dplyr::summarise(
    n_reps = dplyr::n(),
    n_with_divergent = sum(n_divergent > 0),
    pct_with_divergent = round(n_with_divergent / n_reps * 100, 1),
    max_divergent = max(n_divergent),
    .groups = "drop"
  )

cat("\n=== Divergent transitions ===\n")
print(div_report)

# ==========================================================================
# 5. Pass/fail check
# ==========================================================================
cat("\n=== Phase 2 pass/fail ===\n")

for (s_name in unique(results_df$scenario)) {
  s_data <- summary_metrics |> dplyr::filter(scenario == s_name)
  cat(sprintf("\nScenario %s:\n", s_name))

  # Check 1: bias < 0.1
  bias_ok <- abs(s_data$bias) < 0.1
  cat(sprintf("  Bias |%.3f| < 0.1 : %s\n",
              s_data$bias, ifelse(bias_ok, "PASS", "FAIL")))

  # Check 2: coverage 0.85-1.0
  cov_ok <- s_data$coverage_95 >= 0.85
  cat(sprintf("  Coverage %.2f >= 0.85 : %s\n",
              s_data$coverage_95, ifelse(cov_ok, "PASS", "FAIL")))

  # Check 3: divergent < 5%
  div_ok <- s_data$pct_divergent < 5
  cat(sprintf("  Divergent rate %.1f%% < 5%% : %s\n",
              s_data$pct_divergent, ifelse(div_ok, "PASS", "FAIL")))
}

cat("\nIf all PASS, proceed to Phase 3 (real data application).\n")
cat("Output files:\n")
cat("  - outputs/03_summary_metrics.csv\n")
cat("  - outputs/03_recovery_plot.png\n")
