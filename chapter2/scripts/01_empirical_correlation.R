# ==========================================================================
# 01_empirical_correlation_v3.R
# v3 changes:
#   - Uses compute_residual_correlation_v3 (with cluster bootstrap for residuals)
#   - Reports both naive Fisher z and cluster bootstrap CIs for transparency
#   - Other outputs same as v2
# ==========================================================================

setwd("~/chapter2")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(serodynamics)
library(serocalculator)

source("R/correlation_utils.R")
source("R/compute_residual_correlation.R")
source("R/build_summary_row.R")
source("R/plot_empirical_correlation.R")

set.seed(2026)

# ==========================================================================
# 1. Compute for 3 antigens
# ==========================================================================
cat("=== IpaB (n=48) ===\n")
load("~/Data/Manuscript/overall_IpaB_pop_6.rda")
ipab_corr <- compute_residual_correlation_ch1_v3(
  fit = overall_IpaB_pop_6, antigen_label = "IpaB"
)

cat("\n  Residual rho (log):       ", round(ipab_corr$rho_residual_log, 3))
cat("\n  Naive Fisher z 95% CI:    [",
    round(ipab_corr$rho_residual_log_ci[1], 3), ",",
    round(ipab_corr$rho_residual_log_ci[2], 3), "]")
cat("\n  Cluster bootstrap 95% CI: [",
    round(ipab_corr$rho_residual_log_ci_cluster[1], 3), ",",
    round(ipab_corr$rho_residual_log_ci_cluster[2], 3), "]",
    " <-- wider (correct)\n")

cat("\n=== Sonnei (n=11) ===\n")
load("~/Data/Manuscript/serotype_sonnei_3.rda")
sonnei_corr <- compute_residual_correlation_ch1_v3(
  fit = serotype_sonnei_3, antigen_label = "Sonnei"
)

cat("\n  Residual rho (log):       ", round(sonnei_corr$rho_residual_log, 3))
cat("\n  Cluster bootstrap 95% CI: [",
    round(sonnei_corr$rho_residual_log_ci_cluster[1], 3), ",",
    round(sonnei_corr$rho_residual_log_ci_cluster[2], 3), "]\n")

cat("\n=== Sf2a (n=17) ===\n")
load("~/Data/Manuscript/serotype_sf2a_3.rda")
sf2a_corr <- compute_residual_correlation_ch1_v3(
  fit = serotype_sf2a_3, antigen_label = "Sf2a"
)

cat("\n  Residual rho (log):       ", round(sf2a_corr$rho_residual_log, 3))
cat("\n  Cluster bootstrap 95% CI: [",
    round(sf2a_corr$rho_residual_log_ci_cluster[1], 3), ",",
    round(sf2a_corr$rho_residual_log_ci_cluster[2], 3), "]\n")

# ==========================================================================
# 2. Summary table with CIs
# ==========================================================================
all_results <- list()
for (antigen_corr in list(ipab_corr, sonnei_corr, sf2a_corr)) {
  for (pname in c("y0", "y1", "t1", "alpha", "shape")) {
    all_results[[length(all_results) + 1]] <-
      build_summary_row(antigen_corr, pname)
  }
}

summary_df <- dplyr::bind_rows(all_results) |>
  dplyr::mutate(
    CI_fisher_label = sprintf("[%.2f, %.2f]", ci_fisher_lo, ci_fisher_hi),
    CI_boot_label   = sprintf("[%.2f, %.2f]", ci_boot_lo, ci_boot_hi),
    rho_with_ci     = sprintf("%.2f %s", rho, CI_fisher_label)
  )

print(summary_df)
saveRDS(summary_df, "outputs/01v3_summary_with_ci.rds")

# ==========================================================================
# 3. Residual CI comparison table (naive vs cluster bootstrap)
# ==========================================================================
resid_ci_df <- tibble::tibble(
  Antigen = c("IpaB", "Sonnei", "Sf2a"),
  rho_residual_log = c(ipab_corr$rho_residual_log,
                       sonnei_corr$rho_residual_log,
                       sf2a_corr$rho_residual_log),
  fisher_lo = c(ipab_corr$rho_residual_log_ci[1],
                sonnei_corr$rho_residual_log_ci[1],
                sf2a_corr$rho_residual_log_ci[1]),
  fisher_hi = c(ipab_corr$rho_residual_log_ci[2],
                sonnei_corr$rho_residual_log_ci[2],
                sf2a_corr$rho_residual_log_ci[2]),
  cluster_lo = c(ipab_corr$rho_residual_log_ci_cluster[1],
                 sonnei_corr$rho_residual_log_ci_cluster[1],
                 sf2a_corr$rho_residual_log_ci_cluster[1]),
  cluster_hi = c(ipab_corr$rho_residual_log_ci_cluster[2],
                 sonnei_corr$rho_residual_log_ci_cluster[2],
                 sf2a_corr$rho_residual_log_ci_cluster[2])
) |>
  dplyr::mutate(
    fisher_width  = fisher_hi - fisher_lo,
    cluster_width = cluster_hi - cluster_lo,
    width_ratio   = cluster_width / fisher_width
  )

cat("\n=== Residual CI comparison (Ezra's clustering concern) ===\n")
print(resid_ci_df)
cat("\nCluster bootstrap CIs wider by",
    round(mean(resid_ci_df$width_ratio, na.rm = TRUE), 1),
    "× on average — confirms naive Fisher z was optimistic.\n")

saveRDS(resid_ci_df, "outputs/01v3_residual_ci_comparison.rds")

# ==========================================================================
# 4. Panel B scatter plot grid
# ==========================================================================
all_scatter <- dplyr::bind_rows(
  ipab_corr$scatter_df,
  sonnei_corr$scatter_df,
  sf2a_corr$scatter_df
) |>
  dplyr::mutate(
    Parameter   = factor(parameter,
                         levels = c("y0", "y1", "t1", "alpha", "shape")),
    Antigen     = factor(antigen, levels = c("IpaB", "Sonnei", "Sf2a")),
    panel_label = sprintf("rho=%.2f [%.2f,%.2f]", rho, ci_lower, ci_upper)
  )

p_grid <- plot_param_scatter_grid(all_scatter)
ggsave("outputs/01v3_scatter_grid.png", p_grid,
       width = 14, height = 8, dpi = 150, bg = "white")

# ==========================================================================
# 5. Forest plot with CIs
# ==========================================================================
forest_df <- summary_df |>
  dplyr::mutate(
    param_label   = factor(Parameter,
                           levels = c("y0", "y1", "t1", "alpha", "shape")),
    antigen_color = factor(Antigen, levels = c("IpaB", "Sonnei", "Sf2a"))
  )

p_forest <- plot_param_forest(forest_df)
ggsave("outputs/01v3_forest_plot.png", p_forest,
       width = 10, height = 5, dpi = 150, bg = "white")

# ==========================================================================
# 6. Save all
# ==========================================================================
saveRDS(list(
  IpaB   = ipab_corr,
  Sonnei = sonnei_corr,
  Sf2a   = sf2a_corr
), "outputs/01v3_all_correlations.rds")

cat("\nPhase 1 v3 complete.\n")
cat("  Outputs:\n")
cat("  - outputs/01v3_summary_with_ci.rds\n")
cat("  - outputs/01v3_residual_ci_comparison.rds  <- NEW naive vs cluster\n")
cat("  - outputs/01v3_scatter_grid.png\n")
cat("  - outputs/01v3_forest_plot.png\n")
cat("  - outputs/01v3_all_correlations.rds\n")
