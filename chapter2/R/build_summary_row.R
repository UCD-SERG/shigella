#' Build one row of the parameter correlation summary table
#' @param corr output of compute_residual_correlation_ch1_v3()
#' @param param_name one of "y0", "y1", "t1", "alpha", "shape"
#' @return tibble row or NULL
build_summary_row <- function(corr, param_name) {
  r <- corr$param_results[[param_name]]
  if (is.null(r) || is.na(r$rho)) return(NULL)
  tibble::tibble(
    Antigen       = corr$antigen,
    Parameter     = param_name,
    n             = r$n,
    rho           = r$rho,
    ci_fisher_lo  = r$ci_fisher[1],
    ci_fisher_hi  = r$ci_fisher[2],
    ci_boot_lo    = r$ci_boot[1],
    ci_boot_hi    = r$ci_boot[2]
  )
}
