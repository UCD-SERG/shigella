#' Run Phase 0 interactive SLURM reproducibility diagnostic
#'
#' Simulates correlated case data, fits model_2 under an interactive SLURM
#' allocation (`salloc`), extracts diagnostics, and saves a result bundle.
#' Use this to establish a Phase 0 baseline for comparing with Phase 1
#' sbatch results and confirming determinism across allocation modes.
#'
#' @param n Number of subjects to simulate.
#' @param iter_warmup Number of Stan warmup iterations per chain.
#' @param iter_sampling Number of Stan sampling iterations per chain.
#' @param tag File-naming tag, e.g. `"n5"` or `"n48"`.
#' @param output_dir Directory for output files (default `"outputs/phase0"`).
#' @param true_rho_B True Kronecker biomarker correlation (default `0.6`).
#' @param seed Random seed for simulation and Stan (default `20260513`).
#' @param chains Number of MCMC chains (default `2`).
#' @param adapt_delta Stan `adapt_delta` (default `0.95`).
#' @param max_treedepth Stan `max_treedepth` (default `12`).
#' @param compile_dir Directory for compiled Stan binaries. If `NULL`,
#'   defaults to `/tmp/<USER>/cmdstan_bin_phase0_<tag>`.
#' @return Invisibly returns the result bundle list, or `NULL` if the fit
#'   crashed.
#' @export
run_phase0_diagnostic <- function(n,
                                   iter_warmup,
                                   iter_sampling,
                                   tag,
                                   output_dir    = "outputs/phase0",
                                   true_rho_B    = 0.6,
                                   seed          = 20260513L,
                                   chains        = 2L,
                                   adapt_delta   = 0.95,
                                   max_treedepth = 12L,
                                   compile_dir   = NULL) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf(" PHASE 0: INTERACTIVE SLURM REPRODUCIBILITY TEST (%s)\n", tag))
  cat(" Purpose: fit via salloc to compare determinism with Phase 1 sbatch\n")
  cat(strrep("=", 70), "\n\n", sep = "")
  cat(sprintf("Started at:  %s\n", format(Sys.time())))
  cat(sprintf("Host:        %s\n\n", Sys.info()[["nodename"]]))

  dir.create("logs/phase0", recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir,    recursive = TRUE, showWarnings = FALSE)

  status_file <- file.path(output_dir, "PHASE0_STATUS.txt")
  unlink(status_file)
  write_status(status_file, "INIT", "Phase 0 started")

  # ----- 1. Package versions -----
  write_status(status_file, "LOAD_PACKAGES", "logging")
  pkg_versions <- c(
    R            = R.version.string,
    cmdstanr     = as.character(utils::packageVersion("cmdstanr")),
    posterior    = as.character(utils::packageVersion("posterior")),
    shigella     = as.character(utils::packageVersion("shigella")),
    serodynamics = as.character(utils::packageVersion("serodynamics"))
  )
  cmdstan_ver <- tryCatch(
    cmdstanr::cmdstan_version(), error = function(e) "UNKNOWN"
  )
  for (pkg in names(pkg_versions)) {
    cat(sprintf("  %-15s %s\n", pkg, pkg_versions[pkg]))
  }
  cat(sprintf("  %-15s %s\n\n", "cmdstan", cmdstan_ver))
  saveRDS(
    c(pkg_versions, cmdstan = cmdstan_ver),
    file.path(output_dir, "env_versions.rds")
  )
  write_status(status_file, "LOAD_PACKAGES", "OK")

  # ----- 2. Compile dir -----
  write_status(status_file, "COMPILE_DIR", "setting up")
  if (is.null(compile_dir)) {
    user        <- Sys.getenv("USER", unset = "unknown")
    compile_dir <- file.path("/tmp", user,
                             paste0("cmdstan_bin_phase0_", tag))
  }
  if (!dir.exists(compile_dir)) {
    dir.create(compile_dir, recursive = TRUE, mode = "0755")
  }
  cat(sprintf("  compile_dir:    %s\n", compile_dir))
  cat(sprintf("  existing files: %d\n\n", length(list.files(compile_dir))))
  write_status(status_file, "COMPILE_DIR", "OK")

  # ----- 3. Simulate -----
  write_status(status_file, "SIMULATE", "running")
  set.seed(seed)
  omega_B_true <- matrix(c(1, true_rho_B, true_rho_B, 1), 2, 2)
  sim_data <- sim_correlated_case_data(
    n                 = n,
    omega_B           = omega_B_true,
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )
  cat(sprintf("  n_subjects: %d, rows: %d, true rho_B: %.3f\n",
              length(unique(sim_data$id)), nrow(sim_data), true_rho_B))
  saveRDS(sim_data, file.path(output_dir, sprintf("sim_data_%s.rds", tag)))
  write_status(status_file, "SIMULATE", "OK")
  cat("\n")

  # ----- 4. Fit -----
  write_status(status_file, "FIT", "running")
  out_file <- file.path(output_dir, sprintf("one_fit_%s.rds", tag))
  t_start  <- Sys.time()
  saveRDS(
    list(scenario = "phase0_interactive", status = "FIT_STARTED",
         started_at = format(t_start), true_rho_B = true_rho_B, n = n),
    out_file
  )
  fit <- tryCatch({
    run_mod_stan(
      data            = sim_data,
      model           = "model_2",
      chains          = chains,
      iter_warmup     = iter_warmup,
      iter_sampling   = iter_sampling,
      parallel_chains = chains,
      adapt_delta     = adapt_delta,
      max_treedepth   = max_treedepth,
      init            = 0.1,
      with_post       = TRUE,
      compile_dir     = compile_dir,
      refresh         = 100,
      show_messages   = TRUE
    )
  }, error = function(e) {
    cat("\n  [FIT ERROR]:", conditionMessage(e), "\n")
    write_status(status_file, "FIT", paste("CRASHED:", conditionMessage(e)))
    saveRDS(
      list(scenario = "phase0_interactive", status = "FIT_FAILED",
           error = conditionMessage(e), crashed_at = format(Sys.time()),
           true_rho_B = true_rho_B, n = n),
      out_file
    )
    NULL
  })

  elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
  cat(sprintf("\n  Fit elapsed: %.2f min\n\n", elapsed))

  if (is.null(fit)) {
    cat(strrep("!", 70), "\n")
    cat(" PHASE 0 RESULT: FIT CRASHED\n")
    write_status(status_file, "DONE", "Phase 0 FAILED")
    return(invisible(NULL))
  }
  write_status(status_file, "FIT", "OK")

  # ----- 5. Diagnostics -----
  write_status(status_file, "DIAG", "extracting")
  sf          <- attr(fit, "stan_fit")[[1]]
  diag        <- sf$diagnostic_summary(
    diagnostics = c("divergences", "treedepth", "ebfmi")
  )
  total_iters <- chains * iter_sampling
  draws_summary <- tryCatch(
    posterior::summarise_draws(
      sf$draws(variables = "Omega_B[1,2]"),
      "median", "mean", "sd",
      ~ quantile(.x, c(0.025, 0.975), na.rm = TRUE),
      "ess_bulk", "rhat"
    ),
    error = function(e) NULL
  )
  if (!is.null(draws_summary)) {
    cat("\n  Omega_B[1,2] posterior summary:\n")
    print(draws_summary)
  }
  write_status(status_file, "DIAG", "OK")
  cat("\n")

  # ----- 6. Save bundle -----
  write_status(status_file, "SAVE", "writing rds")
  result_bundle <- list(
    scenario           = "phase0_interactive",
    status             = "OK",
    elapsed_min        = elapsed,
    started_at         = format(t_start),
    completed_at       = format(Sys.time()),
    host               = Sys.info()[["nodename"]],
    pkg_versions       = pkg_versions,
    cmdstan_version    = cmdstan_ver,
    true_rho_B         = true_rho_B,
    n_subjects         = n,
    fit_settings       = list(
      chains = chains, warmup = iter_warmup, sampling = iter_sampling,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth
    ),
    diagnostic_summary = diag,
    omega_B_summary    = draws_summary,
    rho_B_posterior    = tryCatch(
      as.vector(posterior::as_draws_array(sf$draws("Omega_B[1,2]"))),
      error = function(e) NULL
    )
  )
  saveRDS(result_bundle, out_file)
  saveRDS(
    result_bundle$diagnostic_summary,
    file.path(output_dir, sprintf("one_fit_%s_diag.rds", tag))
  )
  write_status(status_file, "SAVE", "OK")
  cat(sprintf("  saved -> %s\n\n", out_file))

  # ----- 7. Summary -----
  cat(strrep("=", 70), "\n")
  cat(" PHASE 0 RESULT SUMMARY\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("  Status:              OK\n"))
  cat(sprintf("  Elapsed:             %.2f min\n", elapsed))
  cat(sprintf("  True rho_B:          %+.3f\n",    true_rho_B))
  if (!is.null(draws_summary)) {
    cat(sprintf("  Recovered median:    %+.3f  [%.3f, %.3f]\n",
                draws_summary$median,
                draws_summary$`2.5%`,
                draws_summary$`97.5%`))
    cat(sprintf("  ESS_bulk:            %.0f\n", draws_summary$ess_bulk))
    cat(sprintf("  R-hat:               %.3f\n", draws_summary$rhat))
  }
  cat(sprintf("  Divergent:           %d / %d\n",
              sum(diag$num_divergent), total_iters))
  cat(sprintf("  Max-treedepth hits:  %d / %d\n",
              sum(diag$num_max_treedepth), total_iters))
  cat(strrep("=", 70), "\n\n")

  cat(" NEXT STEP:\n")
  cat(sprintf("   1. Inspect %s + logs/phase0/*.log\n", out_file))
  cat("   2. If divergent rate <= 5% AND R-hat <= 1.01:\n")
  cat("        -> Proceed to Phase 1 (sbatch slurm/phase1_single.sbatch)\n")
  cat("   3. If divergent rate > 10% OR R-hat > 1.02:\n")
  cat("        -> Skip Phase 1-3, jump to Phase 4 diagnosis.\n\n")

  write_status(status_file, "DONE", "Phase 0 completed successfully")
  invisible(result_bundle)
}
