#' Run one simulation replicate for Phase 2
#'
#' Simulates data, fits Stan model_2, extracts rho_B posterior, and returns
#' a result list. Handles errors at each stage with informative status codes.
#'
#' @param s_name character scenario name (e.g. "A", "B", "C")
#' @param scn list with fields: n, rho_B, label
#' @param rep integer replicate index
#' @param n_chains integer number of MCMC chains
#' @param n_iter_warmup integer warmup iterations per chain
#' @param n_iter_sample integer sampling iterations per chain
#' @return named list with fields: scenario, rep, status, and (if OK)
#'   true_rho_B, est_rho_B_median, est_rho_B_mean, est_rho_B_lo,
#'   est_rho_B_hi, bias, n_divergent, elapsed_min
run_one_replicate <- function(s_name, scn, rep,
                              n_chains, n_iter_warmup, n_iter_sample) {
  cat(sprintf("\n[Scenario %s rep %d] %s\n",
              s_name, rep, format(Sys.time())))
  t0 <- Sys.time()

  Omega_B_true <- make_omega_2x2(scn$rho_B)

  set.seed(2026 * 100 + rep)

  # ----- Simulate -----
  sim_dat <- tryCatch({
    sim_correlated_case_data(
      n                 = scn$n,
      omega_B           = Omega_B_true,
      antigen_isos      = c("IgG", "IgA"),
      n_obs_per_subject = 5L
    )
  }, error = function(e) {
    cat(sprintf("  SIM ERROR: %s\n", conditionMessage(e))); NULL
  })
  if (is.null(sim_dat)) {
    return(list(scenario = s_name, rep = rep, status = "SIM_FAILED"))
  }

  # ----- Fit -----
  fit <- tryCatch({
    run_mod_stan(
      data            = sim_dat,
      model           = "model_2",
      chains          = n_chains,
      iter_warmup     = n_iter_warmup,
      iter_sampling   = n_iter_sample,
      parallel_chains = n_chains,
      adapt_delta     = 0.99,
      max_treedepth   = 12,
      init            = 0.1,
      with_post       = TRUE,
      stan_dir        = "inst/stan",
      refresh         = 200,
      show_messages   = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  FIT ERROR: %s\n", conditionMessage(e))); NULL
  })
  if (is.null(fit)) {
    return(list(scenario = s_name, rep = rep, status = "FIT_FAILED"))
  }

  # ----- Extract posterior -----
  rho_B_post <- tryCatch({
    sf <- attr(fit, "stan_fit")[[1]]
    omega_B_draws <- posterior::as_draws_df(
      sf$draws(variables = "Omega_B[1,2]")
    )
    omega_B_draws[["Omega_B[1,2]"]]
  }, error = function(e) {
    cat(sprintf("  EXTRACT ERROR: %s\n", conditionMessage(e))); NULL
  })
  if (is.null(rho_B_post)) {
    return(list(scenario = s_name, rep = rep, status = "EXTRACT_FAILED"))
  }

  # ----- Diagnostics -----
  n_divergent <- tryCatch({
    sf <- attr(fit, "stan_fit")[[1]]
    sum(sf$diagnostic_summary()$num_divergent)
  }, error = function(e) NA_integer_)

  elapsed <- as.numeric(Sys.time() - t0, units = "mins")

  cat(sprintf("  rho_B = %.3f [%.3f, %.3f], bias = %+.3f, %d div, %.1f min\n",
              median(rho_B_post),
              quantile(rho_B_post, 0.025, names = FALSE),
              quantile(rho_B_post, 0.975, names = FALSE),
              median(rho_B_post) - scn$rho_B,
              n_divergent,
              elapsed))

  list(
    scenario          = s_name,
    rep               = rep,
    status            = "OK",
    true_rho_B        = scn$rho_B,
    est_rho_B_median  = median(rho_B_post),
    est_rho_B_mean    = mean(rho_B_post),
    est_rho_B_lo      = quantile(rho_B_post, 0.025, names = FALSE),
    est_rho_B_hi      = quantile(rho_B_post, 0.975, names = FALSE),
    bias              = median(rho_B_post) - scn$rho_B,
    n_divergent       = n_divergent,
    elapsed_min       = elapsed
  )
}
