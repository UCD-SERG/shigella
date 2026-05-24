#' Run Phase 1 SLURM single-job reproducibility diagnostic
#'
#' Simulates correlated case data (identical seed and parameters to Phase 0),
#' fits model_2 inside a SLURM sbatch job, extracts diagnostics, and saves a
#' result bundle. Optionally compares output with a Phase 0 baseline to
#' isolate SLURM-vs-code attribution of any sampling issues.
#'
#' @param n Number of subjects to simulate (must match Phase 0).
#' @param iter_warmup Number of Stan warmup iterations per chain.
#' @param iter_sampling Number of Stan sampling iterations per chain.
#' @param tag File-naming tag, e.g. `"n5"` or `"n48"`.
#' @param output_dir Directory for output files (default `"outputs/phase1"`).
#' @param phase0_dir Directory containing the Phase 0 result bundle
#'   (default `"outputs/phase0"`). Used for the Phase 0 vs Phase 1
#'   comparison table.
#' @param true_rho_B True Kronecker biomarker correlation (default `0.6`).
#' @param seed Random seed - must match Phase 0 (default `20260513`).
#' @param chains Number of MCMC chains (default `2`).
#' @param adapt_delta Stan `adapt_delta` (default `0.95`).
#' @param max_treedepth Stan `max_treedepth` (default `12`).
#' @param compile_dir Directory for compiled Stan binaries. If `NULL`,
#'   falls back to the `STAN_COMPILE_DIR` environment variable, then
#'   `/tmp/<USER>/cmdstan_bin_phase1_<tag>_<jobid>`.
#' @return Invisibly returns the result bundle list, or `NULL` if the fit
#'   crashed.
#' @example inst/examples/run_phase1_diagnostic-examples.R
#' @keywords internal
run_phase1_diagnostic <- function(n,
                                  iter_warmup,
                                  iter_sampling,
                                  tag,
                                  output_dir    = "outputs/phase1",
                                  phase0_dir    = "outputs/phase0",
                                  true_rho_B    = 0.6,
                                  seed          = 20260513L,
                                  chains        = 2L,
                                  adapt_delta   = 0.95,
                                  max_treedepth = 12L,
                                  compile_dir   = NULL) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf(" PHASE 1: SLURM SINGLE JOB DIAGNOSTIC (%s)\n", tag))
  cat(" Purpose: fit inside Slurm; compare with Phase 0 to isolate 
      attribution\n")
  cat(strrep("=", 70), "\n\n", sep = "")

  # ----- 0. SLURM environment -----
  slurm_env <- c(
    SLURM_JOB_ID        = Sys.getenv("SLURM_JOB_ID"),
    SLURM_JOB_NAME      = Sys.getenv("SLURM_JOB_NAME"),
    SLURM_NODELIST      = Sys.getenv("SLURM_NODELIST"),
    SLURM_CPUS_PER_TASK = Sys.getenv("SLURM_CPUS_PER_TASK"),
    SLURM_MEM_PER_NODE  = Sys.getenv("SLURM_MEM_PER_NODE"),
    SLURM_SUBMIT_DIR    = Sys.getenv("SLURM_SUBMIT_DIR"),
    USER                = Sys.getenv("USER"),
    HOSTNAME            = Sys.info()[["nodename"]],
    TMPDIR              = Sys.getenv("TMPDIR")
  )
  cat("=== SLURM environment ===\n")
  for (env_name in names(slurm_env)) {
    cat(sprintf("  %-22s = %s\n", env_name, slurm_env[env_name]))
  }
  cat("\n")
  cat(sprintf("Started at:  %s\n", format(Sys.time())))
  cat(sprintf("R version:   %s\n\n", R.version.string))

  dir.create("logs/phase1", recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir,    recursive = TRUE, showWarnings = FALSE)

  job_id <- slurm_env["SLURM_JOB_ID"]
  if (job_id == "") job_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  status_file <- file.path(output_dir,
                           sprintf("PHASE1_STATUS_%s.txt", job_id))
  write_status(status_file, "INIT",
               sprintf("Phase 1 started, jobid=%s", job_id))

  # ----- 1. Package versions -----
  write_status(status_file, "LOAD_PACKAGES", "logging")
  pkg_versions <- c(
    R            = R.version.string,
    cmdstanr     = .safe_pkg_version("cmdstanr"),
    posterior    = .safe_pkg_version("posterior"),
    shigella     = .safe_pkg_version("shigella"),
    serodynamics = .safe_pkg_version("serodynamics")
  )
  cmdstan_ver <- tryCatch(
    cmdstanr::cmdstan_version(), error = function(e) "UNKNOWN"
  )
  cat("=== Package versions ===\n")
  for (pkg in names(pkg_versions)) {
    cat(sprintf("  %-15s %s\n", pkg, pkg_versions[pkg]))
  }
  cat(sprintf("  %-15s %s\n\n", "cmdstan", cmdstan_ver))
  write_status(status_file, "LOAD_PACKAGES", "OK")

  # ----- 2. Compile dir (SLURM-specific, per-task subdir) -----
  write_status(status_file, "COMPILE_DIR", "setting up")
  if (is.null(compile_dir)) {
    compile_dir <- Sys.getenv("STAN_COMPILE_DIR", unset = "")
    if (compile_dir == "") {
      user        <- Sys.getenv("USER", unset = "unknown")
      compile_dir <- file.path(
        "/tmp", user, sprintf("cmdstan_bin_phase1_%s_%s", tag, job_id)
      )
    }
  }
  if (!dir.exists(compile_dir)) {
    dir.create(compile_dir, recursive = TRUE, mode = "0755")
  }
  cat(sprintf("=== compile_dir: %s ===\n",    compile_dir))
  cat(sprintf("    existing files: %d\n\n", length(list.files(compile_dir))))
  write_status(status_file, "COMPILE_DIR", "OK")

  # ----- 3. Simulate (identical to Phase 0 - same seed, same n) -----
  write_status(status_file, "SIMULATE", "running")
  set.seed(seed)
  omega_B_true <- matrix(c(1, true_rho_B, true_rho_B, 1), 2, 2)
  sim_data <- sim_correlated_case_data(
    n                 = n,
    omega_B           = omega_B_true,
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )
  cat(sprintf("=== Simulated: n=%d, rho_B=%.1f, total rows=%d ===\n\n",
              length(unique(sim_data$id)), true_rho_B, nrow(sim_data)))
  write_status(status_file, "SIMULATE", "OK")

  # ----- 4. Save STARTED placeholder -----
  scenario <- sprintf("phase1_slurm_single_%s", tag)
  out_file <- file.path(output_dir,
                        sprintf("one_fit_%s_jobid_%s.rds", tag, job_id))
  saveRDS(
    list(scenario = scenario, status = "FIT_STARTED", job_id = job_id,
         started_at = format(Sys.time()), slurm_env = slurm_env,
         pkg_versions = pkg_versions, cmdstan_version = cmdstan_ver,
         true_rho_B = true_rho_B, n = n),
    out_file
  )

  # ----- 5. Fit -----
  write_status(status_file, "FIT", "running")
  cat("=== Fitting model_2 (Kronecker) ===\n")
  t_start <- Sys.time()
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
      list(scenario = scenario, status = "FIT_FAILED", job_id = job_id,
           error = conditionMessage(e), crashed_at = format(Sys.time()),
           slurm_env = slurm_env, pkg_versions = pkg_versions),
      out_file
    )
    NULL
  })

  elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
  cat(sprintf("\n=== Fit elapsed: %.2f min ===\n\n", elapsed))

  if (is.null(fit)) {
    phase0_file <- file.path(phase0_dir, sprintf("one_fit_%s.rds", tag))
    cat(strrep("!", 70), "\n")
    cat(" PHASE 1 RESULT: FIT CRASHED INSIDE SLURM\n")
    cat(sprintf(" Compare with %s to determine:\n", phase0_file))
    cat("   - If Phase 0 OK but Phase 1 FAIL -> Slurm env issue\n")
    cat("   - If both fail -> code / model identifiability issue\n")
    cat(strrep("!", 70), "\n")
    write_status(status_file, "DONE", "Phase 1 FAILED")
    return(invisible(NULL))
  }
  write_status(status_file, "FIT", "OK")

  # ----- 6. Diagnostics -----
  write_status(status_file, "DIAG", "extracting")
  sf          <- attr(fit, "stan_fit")[[1]]
  diag        <- sf$diagnostic_summary(
    diagnostics = c("divergences", "treedepth", "ebfmi")
  )
  total_iters  <- chains * iter_sampling
  draws_summary <- tryCatch(
    posterior::summarise_draws(
      sf$draws(variables = "Omega_B[1,2]"),
      "median", "mean", "sd",
      ~ quantile(.x, c(0.025, 0.975), na.rm = TRUE),
      "ess_bulk", "rhat"
    ),
    error = function(e) NULL
  )
  cat("=== Diagnostics ===\n")
  cat(sprintf("  divergent:     %d / %d (%.2f%%)\n",
              sum(diag$num_divergent), total_iters,
              100 * sum(diag$num_divergent) / total_iters))
  cat(sprintf("  max-treedepth: %d / %d (%.2f%%)\n",
              sum(diag$num_max_treedepth), total_iters,
              100 * sum(diag$num_max_treedepth) / total_iters))
  cat(sprintf("  E-BFMI:        %s\n",
              paste(sprintf("%.3f", diag$ebfmi), collapse = ", ")))
  if (!is.null(draws_summary)) {
    cat("\n  Omega_B[1,2] summary:\n")
    print(draws_summary)
  } else {
    cat("\n  [INFO] Omega_B[1,2] not available in this fit.\n")
  }
  write_status(status_file, "DIAG", "OK")

  # ----- 7. Save bundle -----
  write_status(status_file, "SAVE", "writing rds")
  result_bundle <- list(
    scenario           = scenario,
    status             = "OK",
    job_id             = job_id,
    elapsed_min        = elapsed,
    started_at         = format(t_start),
    completed_at       = format(Sys.time()),
    slurm_env          = slurm_env,
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

  # ----- 8. Compare with Phase 0 if available -----
  phase0_file <- file.path(phase0_dir, sprintf("one_fit_%s.rds", tag))
  if (file.exists(phase0_file)) {
    ph0 <- readRDS(phase0_file)
    if (!is.null(ph0$omega_B_summary) && !is.null(draws_summary)) {
      cat("\n=== Phase 0 vs Phase 1 comparison ===\n")
      cmp <- data.frame(
        metric = c("status", "elapsed_min", "post_median",
                   "post_lo_2.5", "post_hi_97.5",
                   "ess_bulk", "rhat", "n_divergent", "n_treedepth"),
        phase0 = c(
          ph0$status,
          round(ph0$elapsed_min, 2),
          round(ph0$omega_B_summary$median, 3),
          round(ph0$omega_B_summary$`2.5%`, 3),
          round(ph0$omega_B_summary$`97.5%`, 3),
          round(ph0$omega_B_summary$ess_bulk, 0),
          round(ph0$omega_B_summary$rhat, 3),
          sum(ph0$diagnostic_summary$num_divergent),
          sum(ph0$diagnostic_summary$num_max_treedepth)
        ),
        phase1 = c(
          "OK",
          round(elapsed, 2),
          round(draws_summary$median, 3),
          round(draws_summary$`2.5%`, 3),
          round(draws_summary$`97.5%`, 3),
          round(draws_summary$ess_bulk, 0),
          round(draws_summary$rhat, 3),
          sum(diag$num_divergent),
          sum(diag$num_max_treedepth)
        )
      )
      print(cmp, row.names = FALSE)
      saveRDS(cmp, file.path(output_dir,
                             sprintf("p0_vs_p1_comparison_%s.rds", job_id)))
    }
  } else {
    cat("\n  [INFO] Phase 0 result not found - comparison skipped.\n")
    cat(sprintf("         Run phase0_interactive_reproducibility_%s.R first",
                tag))
    cat(" for direct comparison.\n")
  }

  write_status(status_file, "SAVE", "OK")
  cat(sprintf("\n=== Phase 1 complete ===\n    Results: %s\n\n", out_file))
  write_status(status_file, "DONE", "Phase 1 OK")
  invisible(result_bundle)
}
