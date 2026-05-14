#' @title Prepare data for Stan backend
#' @description
#' Converts the output of [serodynamics::prep_data()] (a
#' `prepped_jags_data` list) into the list format required by the Stan
#' models `model_1.stan` and `model_2.stan`. Handles NA-padding of the
#' ragged observation array.
#'
#' The Stan models expect:
#' - `N`: number of subjects
#' - `K`: number of antigen-isotype biomarkers
#' - `P`: number of kinetic parameters
#' - `max_obs`: max number of observations per subject
#' - `n_obs[N]`: actual number of observations per subject
#'   (ragged handling)
#' - `time_obs[N, max_obs]`: observation times (NA -> 0, ignored by
#'   likelihood)
#' - `log_y[N, max_obs, K]`: log-transformed antibody observations
#'
#' @param prepped_jags_data output from
#'   [serodynamics::prep_data()] — a list with elements `smpl.t`,
#'   `logy`, `nsmpl`, `nsubj`, `n_antigen_isos`.
#' @param drop_newperson [logical] whether to drop the JAGS dummy
#'   "newperson" row before passing to Stan. Default `TRUE` because Stan
#'   handles posterior prediction through the `generated quantities`
#'   block rather than through a missing-data dummy subject.
#'
#' @returns a named [list] ready to pass to a compiled Stan model via
#'   `cmdstanr::cmdstan_model()$sample()`.
#' @export
#' @example inst/examples/prep_data_stan-examples.R
prep_data_stan <- function(prepped_jags_data,
                           drop_newperson = TRUE) {

  if (!inherits(prepped_jags_data, "prepped_jags_data")) {
    cli::cli_abort(c(
      "{.arg prepped_jags_data} must be a {.cls prepped_jags_data} object",
      "i" = "Did you forget to call {.fn prep_data} first?"
    ))
  }

  # Extract arrays
  smpl_t  <- prepped_jags_data$smpl.t    # [nsubj, max_visits]
  logy    <- prepped_jags_data$logy      # [nsubj, max_visits, K]
  nsmpl   <- as.integer(prepped_jags_data$nsmpl)
  K       <- prepped_jags_data$n_antigen_isos
  N_full  <- prepped_jags_data$nsubj

  ids_all <- attr(prepped_jags_data, "ids")

  # Drop the "newperson" dummy row added by prep_data()
  if (drop_newperson && "newperson" %in% ids_all) {
    keep_idx <- which(ids_all != "newperson")
    smpl_t   <- smpl_t[keep_idx, , drop = FALSE]
    logy     <- logy[keep_idx, , , drop = FALSE]
    nsmpl    <- nsmpl[keep_idx]
    ids_kept <- ids_all[keep_idx]
    N        <- length(keep_idx)
  } else {
    ids_kept <- ids_all
    N        <- N_full
  }

  max_obs <- ncol(smpl_t)
  P       <- 5L   

  # Replace NA with 0 (time) and 0 (log_y) — Stan ignores these via n_obs[i] guard
  # in the likelihood loop (for (t_idx in 1:n_obs[i])).
  time_obs <- smpl_t
  time_obs[is.na(time_obs)] <- 0

  log_y <- logy
  log_y[is.na(log_y)] <- 0

  # Sanity checks
  if (any(nsmpl > max_obs)) {
    cli::cli_abort(
      "n_obs[{which(nsmpl > max_obs)}] > max_obs. Array sizes inconsistent."
    )
  }
  if (any(nsmpl == 0)) {
    cli::cli_warn(
      "Subject(s) with 0 observations detected; these contribute no likelihood."
    )
  }

  antigens <- attr(prepped_jags_data, "antigens")

  stan_data <- list(
    N        = N,
    K        = as.integer(K),
    P        = P,
    max_obs  = as.integer(max_obs),
    n_obs    = nsmpl,
    time_obs = time_obs,
    log_y    = log_y
  )

  # Attach metadata for postprocessing
  attr(stan_data, "ids")      <- ids_kept
  attr(stan_data, "antigens") <- antigens

  return(stan_data)
}
