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

source("R/summarize_sim_results.R")
source("R/plot_recovery.R")

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
summary_metrics <- summarize_sim_results(results_df)

cat("=== Summary metrics ===\n")
print(summary_metrics)

saveRDS(summary_metrics, "outputs/03_summary_metrics.rds")
write.csv(summary_metrics, "outputs/03_summary_metrics.csv", row.names = FALSE)

# ==========================================================================
# 3. Recovery plot
# ==========================================================================
p_recovery <- plot_recovery(results_df)

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
report_pass_fail(summary_metrics, results_df)

cat("Output files:\n")
cat("  - outputs/03_summary_metrics.csv\n")
cat("  - outputs/03_recovery_plot.png\n")
