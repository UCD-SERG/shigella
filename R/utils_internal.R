# ============================================================
# utils_internal.R
# Internal helper functions that replace external ::: calls.
#
# These functions replicate the minimal subset of serodynamics
# and serocalculator internals needed by the shigella package,
# removing the need for triple-colon (:::) imports.
#
# Source reference:
#   serodynamics::ab()  — two-phase antibody kinetic model
#   serodynamics::bt()  — power-law decay kernel
#   serodynamics:::get_timeindays_var()
#   serodynamics:::use_att_names()
#   serocalculator::get_values_var()
#   serocalculator::as_case_data()
# ============================================================


# ── Two-phase antibody kinetic model ────────────────────────────────────────

#' Power-Law Decay Kernel (Internal)
#'
#' Computes the decay factor for the power-law waning phase of the
#' two-phase antibody kinetic model.
#'
#' @param t Numeric vector of times (days since infection).
#' @param t1 Numeric scalar or vector: time to peak (days).
#' @param alpha Numeric: decay rate (day^{-1}).
#' @param shape Numeric: decay shape parameter (rho).
#'
#' @return Numeric vector of decay factors (multiplied by peak to get
#'   antibody level during decay phase).
#'
#' @keywords internal
#' @noRd
.bt <- function(t, t1, alpha, shape) {
  (1 + (shape - 1) * alpha * (t - t1))^(-1 / (shape - 1))
}


#' Two-Phase Antibody Trajectory (Internal)
#'
#' Evaluates the closed-form two-phase rise-peak-decay antibody model
#' at given times. Replaces `serodynamics:::ab()`.
#'
#' @param t Numeric vector of times (days since symptom onset).
#' @param y0 Log-scale baseline: actual baseline = exp(y0).
#' @param y1 Log-scale (peak - baseline): actual peak = exp(y0) + exp(y1).
#' @param t1 Log-scale time to peak: actual t1 = exp(t1) days.
#' @param alpha Log-scale decay rate: actual alpha = exp(alpha) day^{-1}.
#' @param shape Log-scale (shape - 1): actual shape = exp(shape) + 1.
#'
#' @return Numeric vector of predicted antibody levels (natural MFI scale).
#'
#' @details
#' All five parameters are on the **log scale** as stored by the JAGS model:
#' \itemize{
#'   \item Baseline = exp(y0)
#'   \item Peak     = exp(y0) + exp(y1)
#'   \item t1       = exp(t1)
#'   \item alpha    = exp(alpha)
#'   \item shape    = exp(shape) + 1
#' }
#' The model is:
#' \itemize{
#'   \item Rise phase (t <= exp(t1)): exponential growth from baseline to peak
#'   \item Decay phase (t > exp(t1)): power-law waning governed by alpha and shape
#' }
#'
#' @keywords internal
#' @noRd
.ab <- function(t, y0, y1, t1, alpha, shape) {
  # Transform from log scale to natural scale
  y0_nat    <- exp(y0)
  peak      <- exp(y0) + exp(y1)
  t1_nat    <- exp(t1)
  alpha_nat <- exp(alpha)
  shape_nat <- exp(shape) + 1

  # Rise phase: exponential growth
  mu_y <- log(peak / y0_nat) / t1_nat
  rise <- y0_nat * exp(mu_y * t)

  # Decay phase: power-law waning
  decay <- peak * .bt(t, t1_nat, alpha_nat, shape_nat)

  # Select phase based on time
  ifelse(t <= t1_nat, rise, decay)
}


# ── Case-data attribute helpers ─────────────────────────────────────────────

#' Get the Time-in-Days Variable Name from a Case-Data Object (Internal)
#'
#' Replaces `serodynamics:::get_timeindays_var()`.
#'
#' @param dataset A case_data object.
#' @return Character scalar: the column name storing time in days.
#' @keywords internal
#' @noRd
.get_timeindays_var <- function(dataset) {
  v <- attr(dataset, "timeindays")
  if (!is.null(v)) return(v)
  # Fallback detection
  candidates <- c("timeindays", "Actual day", "timepoint")
  found <- intersect(candidates, names(dataset))
  if (length(found) > 0) return(found[1])
  stop("Cannot determine time variable from dataset.")
}


#' Get the Value Variable Name from a Case-Data Object (Internal)
#'
#' Replaces `serocalculator::get_values_var()`.
#'
#' @param dataset A case_data object.
#' @return Character scalar: the column name storing antibody values.
#' @keywords internal
#' @noRd
.get_values_var <- function(dataset) {
  v <- attr(dataset, "value_var")
  if (!is.null(v)) return(v)
  if ("result" %in% names(dataset)) return("result")
  mfi_cols <- names(dataset)[grepl("_MFI$", names(dataset))]
  if (length(mfi_cols) == 1) return(mfi_cols)
  stop("Cannot determine value variable from dataset.")
}


#' Rename Columns Using Case-Data Attributes (Internal)
#'
#' Replaces `serodynamics:::use_att_names()`. Renames id, time, and value
#' columns to standardised names (`Subject`, `t`, `result`).
#'
#' @param dataset A case_data object.
#' @return The dataset with renamed columns.
#' @keywords internal
#' @noRd
.use_att_names <- function(dataset) {
  id_var    <- attr(dataset, "id_var") %||% "id"
  time_var  <- .get_timeindays_var(dataset)
  value_var <- .get_values_var(dataset)
  iso_var   <- attr(dataset, "biomarker_var") %||% "antigen_iso"

  out <- dataset
  if (id_var != "Subject" && id_var %in% names(out)) {
    out[["Subject"]] <- out[[id_var]]
  }
  if (iso_var != "Iso_type" && iso_var %in% names(out)) {
    out[["Iso_type"]] <- out[[iso_var]]
  }
  if (time_var != "t" && time_var %in% names(out)) {
    out[["t"]] <- out[[time_var]]
  }
  if (value_var != "result" && value_var %in% names(out)) {
    out[["result"]] <- out[[value_var]]
  }
  out
}


#' Create a Case-Data Object (Internal)
#'
#' Lightweight replacement for `serocalculator::as_case_data()`.
#' Attaches the required attributes and sets the class.
#'
#' @param data A data frame.
#' @param id_var Column name for individual IDs.
#' @param biomarker_var Column name for isotype/biomarker labels.
#' @param time_in_days Column name for time in days.
#' @param value_var Column name for antibody values.
#'
#' @return A tibble with class `c("case_data", "tbl_df", "tbl", "data.frame")`
#'   and the four attribute slots set.
#'
#' @keywords internal
#' @noRd
.as_case_data <- function(data,
                          id_var        = "id",
                          biomarker_var = "antigen_iso",
                          time_in_days  = "timeindays",
                          value_var     = "result") {
  out <- tibble::as_tibble(data)
  # Ensure the id column is also available as "id"
  if (id_var != "id" && !"id" %in% names(out)) {
    out[["id"]] <- out[[id_var]]
  }
  class(out) <- c("case_data", class(out))
  attr(out, "id_var")        <- id_var
  attr(out, "biomarker_var") <- biomarker_var
  attr(out, "timeindays")    <- time_in_days
  attr(out, "value_var")     <- value_var
  out
}


#' Null-coalescing operator (Internal)
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
