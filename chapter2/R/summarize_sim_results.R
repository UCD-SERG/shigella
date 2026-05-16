#' Summarize Phase 2 simulation results
#'
#' @param results_df data.frame of successful simulation results
#'   (already filtered to status == "OK")
#' @return tibble of summary metrics per scenario
summarize_sim_results <- function(results_df) {
  results_df |>
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
}

#' Report pass/fail checks for Phase 2 simulation
#'
#' @param summary_metrics tibble returned by summarize_sim_results()
#' @param results_df data.frame of successful simulation results
#' @return invisible NULL (called for side effects: printing)
report_pass_fail <- function(summary_metrics, results_df) {
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
  invisible(NULL)
}
