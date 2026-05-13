#' Chapter 1 fit — residual & parameter correlation (v3)
#'
#'
#' @param fit sr_model object (IgG + IgA both present)
#' @param antigen_label example "IpaB"
#' @param n_boot bootstrap replicates (default 1000)
#'
#' @section Fitted value calculation (Ezra's clarification):
#' Sam's calc_fit_mod uses PLUG-IN estimator:
#'   y_hat = f(t; median_s(theta^s))
#' This differs from posterior-median-of-fitted:
#'   y_alt = median_s(f(t; theta^s))
#' by Jensen's inequality (two-phase curve is nonlinear).
#' Phase 1 diagnostic OK with plug-in; Phase 2 Stan uses full posterior.
#'
#' @section Scale of residuals:
#' Sam stores residuals natural-scale: residual = observed - fitted.
#' We reconstruct log-scale residual for analysis because Ch1 likelihood
#' is on log MFI: logy ~ dnorm(log(mu), tau.logy).
#'
compute_residual_correlation_ch1_v3 <- function(fit,
                                                 antigen_label = "IpaB",
                                                 n_boot = 1000) {

  fr <- attr(fit, "fitted_residuals")
  if (is.null(fr)) {
    stop("fit has no 'fitted_residuals' attribute.")
  }

  fr_igg <- fr |> dplyr::filter(Iso_type == "IgG")
  fr_iga <- fr |> dplyr::filter(Iso_type == "IgA")

  # ===== Residual correlation =====
  merged <- dplyr::inner_join(
    fr_igg |> dplyr::select(Subject, t, residual_igg = residual,
                             fitted_igg = fitted),
    fr_iga |> dplyr::select(Subject, t, residual_iga = residual,
                             fitted_iga = fitted),
    by = c("Subject", "t")
  ) |>
    dplyr::filter(!is.na(residual_igg), !is.na(residual_iga))

  if (nrow(merged) < 5) {
    warning("Only ", nrow(merged), " paired observations")
    return(NULL)
  }

  rho_residual <- cor(merged$residual_igg, merged$residual_iga)

  # Reconstruct log-scale residuals
  # residual = observed - fitted (natural scale)
  # log_resid = log(observed) - log(fitted) = log((fitted+residual)/fitted)
  merged_log <- merged |>
    dplyr::mutate(
      obs_igg = fitted_igg + residual_igg,
      obs_iga = fitted_iga + residual_iga,
      log_resid_igg = log(pmax(obs_igg, 0.01)) - log(pmax(fitted_igg, 0.01)),
      log_resid_iga = log(pmax(obs_iga, 0.01)) - log(pmax(fitted_iga, 0.01))
    ) |>
    dplyr::filter(is.finite(log_resid_igg), is.finite(log_resid_iga))

  rho_residual_log <- cor(merged_log$log_resid_igg, merged_log$log_resid_iga)

  rho_residual_log_ci <- cor.test(rho_residual_log, nrow(merged_log))

  # Cluster bootstrap CI (subject-level resample) — more defensible
  rho_residual_log_ci_cluster <- cluster_bootstrap_residual_ci(
    merged_log, n_boot = n_boot
  )

  # ===== Parameter-level correlation + CI =====
  # NOTE: Parameter correlations are on subject-level (each subject gives one
  # pair of IgG/IgA medians), so observations ARE independent across subjects.
  # Fisher z CI here is valid.

  extract_param_medians <- function(fit_obj, iso_label) {
    fit_obj |>
      dplyr::filter(Iso_type == iso_label,
                    Parameter %in% c("y0", "y1", "t1", "alpha", "shape")) |>
      dplyr::group_by(Subject, Parameter) |>
      dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")
  }

  med_igg <- extract_param_medians(fit, "IgG")
  med_iga <- extract_param_medians(fit, "IgA")

  med_wide <- dplyr::full_join(
    med_igg |> dplyr::rename(IgG = med),
    med_iga |> dplyr::rename(IgA = med),
    by = c("Subject", "Parameter")
  )

  param_results <- list()
  scatter_data <- list()

  for (pname in c("y0", "y1", "t1", "alpha", "shape")) {
    sub <- med_wide |>
      dplyr::filter(Parameter == pname) |>
      dplyr::filter(!is.na(IgG), !is.na(IgA))

    if (nrow(sub) < 5) {
      param_results[[pname]] <- list(rho = NA, ci_fisher = c(NA, NA),
                                      ci_boot = c(NA, NA), n = nrow(sub))
      next
    }

    rho_hat <- cor(sub$IgG, sub$IgA)

    # Fisher z CI 
    ci_fisher <- fisher_z_ci(rho_hat, nrow(sub))

    # Bootstrap CI 
    ci_boot <- bootstrap_cor_ci(sub$IgG, sub$IgA, n_boot = n_boot)

    param_results[[pname]] <- list(
      rho = rho_hat,
      ci_fisher = ci_fisher,
      ci_boot = ci_boot,
      n = nrow(sub)
    )

    scatter_data[[pname]] <- sub |>
      dplyr::mutate(
        antigen = antigen_label,
        parameter = pname,
        rho = rho_hat,
        ci_lower = ci_fisher[1],
        ci_upper = ci_fisher[2]
      )
  }

  scatter_df <- dplyr::bind_rows(scatter_data)

  list(
    antigen                    = antigen_label,
    n_paired_obs               = nrow(merged),
    n_subjects                 = length(unique(merged$Subject)),
    rho_residual               = rho_residual,
    rho_residual_log           = rho_residual_log,
    rho_residual_log_ci        = rho_residual_log_ci,           # optimistic
    rho_residual_log_ci_cluster = rho_residual_log_ci_cluster,  # defensible
    param_results              = param_results,
    scatter_df                 = scatter_df,
    merged_residuals           = merged,
    fitted_value_method = paste0(
      "Plug-in: fitted = f(t; median_s(theta^s)) per Sam's calc_fit_mod. ",
      "Differs from posterior-median-of-fitted by Jensens inequality since ",
      "two-phase curve is nonlinear. OK for Phase 1 diagnostic."
    ),
    ci_method_notes = paste0(
      "Parameter CI: Fisher z valid (subject-level independence). ",
      "Residual CI (naive Fisher z): optimistic due to within-subject clustering. ",
      "Cluster bootstrap CI also provided for defensible reporting."
    )
  )
}


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
    # Resample subjects with replacement
    sampled_subjects <- sample(subjects, n_subj, replace = TRUE)

    # Build resampled dataset
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


#' LR test for residual independence
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
