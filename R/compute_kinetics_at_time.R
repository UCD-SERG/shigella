# Helper: compute log-scale mean for one biomarker at one time point.
# Implements the two-phase power-law kinetics model.
#' @keywords internal
#' @noRd
.compute_kinetics_at_time <- function(
    log_y0, log_y1m0, # nolint: object_name_linter. Stan param name.
    log_t1, log_alpha, log_rm1, tt) {
  y0    <- exp(log_y0)
  y1    <- y0 + exp(log_y1m0)
  t1_j  <- exp(log_t1)
  alpha <- exp(log_alpha)
  shape <- exp(log_rm1) + 1

  if (tt <= t1_j) {
    beta_growth <- (log(y1) - log(y0)) / t1_j
    return(log(y0) + beta_growth * tt)
  }

  # Tolerance is sqrt(.Machine$double.eps) (~1.5e-8): well outside the
  # production prior on shape (typical |1 - shape| > 0.3), but excludes
  # the numerically-unstable region where log(term) / (1 - shape) loses
  # meaningful precision.
  if (abs(1 - shape) < sqrt(.Machine$double.eps)) {
    cli::cli_abort(c(
      "shape ~= 1 is degenerate for the two-phase model",
      "i" = "log_rm1 = {log_rm1} produces shape = {shape}, |1 - shape| = 
      {abs(1 - shape)}",
      "i" = paste0("the decay-phase formula log(term) / (1 - shape) is",
                   " undefined or numerically unstable in this region.")
    ))
  }

  term <- y1^(1 - shape) - (1 - shape) * alpha * (tt - t1_j)
  if (term <= 0) {
    cli::cli_abort(c(
      "Trajectory infeasibility detected at t = {tt}",
      "i" = paste0("term = {term} <= 0; parameter combination",
                   " is invalid for the two-phase model")
    ))
  }
  log(term) / (1 - shape)
}
