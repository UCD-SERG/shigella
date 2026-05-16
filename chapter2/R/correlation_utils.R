# Statistical utilities for correlation analysis (Chapter 2)

#' Fisher z-transformation CI for Pearson correlation
#' @param rho sample correlation
#' @param n sample size (NB: assumes independent observations)
#' @param alpha confidence level
fisher_z_ci <- function(rho, n, alpha = 0.05) {
  if (n < 4 || abs(rho) >= 1) return(c(NA, NA))
  z <- 0.5 * log((1 + rho) / (1 - rho))
  se <- 1 / sqrt(n - 3)
  crit <- qnorm(1 - alpha / 2)
  c(
    lower = (exp(2 * (z - crit * se)) - 1) / (exp(2 * (z - crit * se)) + 1),
    upper = (exp(2 * (z + crit * se)) - 1) / (exp(2 * (z + crit * se)) + 1)
  )
}


#' Bootstrap percentile CI for Pearson correlation (assumes iid)
bootstrap_cor_ci <- function(x, y, n_boot = 1000, alpha = 0.05) {
  n <- length(x)
  if (n < 5) return(c(NA, NA))
  boot_rhos <- replicate(n_boot, {
    idx <- sample(n, n, replace = TRUE)
    xi <- x[idx]; yi <- y[idx]
    if (sd(xi) == 0 || sd(yi) == 0) return(NA)
    cor(xi, yi)
  })
  boot_rhos <- boot_rhos[!is.na(boot_rhos)]
  c(
    lower = quantile(boot_rhos, alpha / 2, names = FALSE),
    upper = quantile(boot_rhos, 1 - alpha / 2, names = FALSE)
  )
}


#' Cluster bootstrap CI for residual correlation
#'
#' Resamples SUBJECTS (not observations) to preserve within-subject clustering.
#' This gives properly calibrated CI for log-scale residual correlation.
#'
#' @param merged_log data with columns Subject, log_resid_igg, log_resid_iga
#' @param n_boot number of bootstrap replicates
#' @param alpha confidence level
cluster_bootstrap_residual_ci <- function(merged_log, n_boot = 1000,
                                          alpha = 0.05) {
  subjects <- unique(merged_log$Subject)
  n_subj <- length(subjects)

  if (n_subj < 5) return(c(NA, NA))

  boot_rhos <- replicate(n_boot, {
    sampled_subjects <- sample(subjects, n_subj, replace = TRUE)

    resampled <- do.call(rbind, lapply(sampled_subjects, function(s) {
      merged_log[merged_log$Subject == s, ]
    }))

    if (nrow(resampled) < 5 ||
        sd(resampled$log_resid_igg) == 0 ||
        sd(resampled$log_resid_iga) == 0) {
      return(NA)
    }
    cor(resampled$log_resid_igg, resampled$log_resid_iga)
  })

  boot_rhos <- boot_rhos[!is.na(boot_rhos)]
  c(
    lower = quantile(boot_rhos, alpha / 2, names = FALSE),
    upper = quantile(boot_rhos, 1 - alpha / 2, names = FALSE)
  )
}


#' Extract posterior median per subject/parameter for one isotype
#'
#' @param fit_obj sr_model object
#' @param iso_label character, e.g. "IgG" or "IgA"
extract_param_medians <- function(fit_obj, iso_label) {
  fit_obj |>
    dplyr::filter(Iso_type == iso_label,
                  Parameter %in% c("y0", "y1", "t1", "alpha", "shape")) |>
    dplyr::group_by(Subject, Parameter) |>
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")
}


#' LR test for residual independence (chi-squared with df=1)
lr_test_independence <- function(rho, n) {
  if (abs(rho) >= 1 || n < 4) return(list(statistic = NA, p_value = NA))
  lambda <- n * log(1 / (1 - rho^2))
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  list(
    statistic = lambda,
    p_value = p_value,
    reject_H0 = lambda > 3.84
  )
}
