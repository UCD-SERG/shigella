# ==========================================================================
# phase1_single_diagnostic.R
# ==========================================================================

setwd("~/shigella/chapter2")

cat("\n", strrep("=", 70), "\n", sep = "")
cat(" PHASE 1: SLURM SINGLE JOB DIAGNOSTIC\n")
cat(" Purpose: Run identical fit INSIDE Slurm (single task)\n")
cat(" Compare with Phase 0 to isolate SLURM-vs-code attribution\n")
cat(strrep("=", 70), "\n\n", sep = "")

# ----- 0. Capture SLURM environment -----
slurm_env <- c(
  SLURM_JOB_ID       = Sys.getenv("SLURM_JOB_ID"),
  SLURM_JOB_NAME     = Sys.getenv("SLURM_JOB_NAME"),
  SLURM_NODELIST     = Sys.getenv("SLURM_NODELIST"),
  SLURM_CPUS_PER_TASK= Sys.getenv("SLURM_CPUS_PER_TASK"),
  SLURM_MEM_PER_NODE = Sys.getenv("SLURM_MEM_PER_NODE"),
  SLURM_SUBMIT_DIR   = Sys.getenv("SLURM_SUBMIT_DIR"),
  USER               = Sys.getenv("USER"),
  HOSTNAME           = Sys.info()[["nodename"]],
  TMPDIR             = Sys.getenv("TMPDIR")
)

cat("=== SLURM environment ===\n")
for (n in names(slurm_env)) cat(sprintf("  %-22s = %s\n", n, slurm_env[n]))
cat("\n")

cat(sprintf("Started at:  %s\n", format(Sys.time())))
cat(sprintf("R version:   %s\n\n", R.version.string))

dir.create("logs/phase1",    recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/phase1", recursive = TRUE, showWarnings = FALSE)

# Status file — uses job id so multiple submissions don't clobber
job_id <- slurm_env["SLURM_JOB_ID"]
if (job_id == "") job_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
status_file <- sprintf("outputs/phase1/PHASE1_STATUS_%s.txt", job_id)

write_status(status_file, "INIT", sprintf("Phase 1 started, jobid=%s", job_id))

# ----- 1. Load packages -----
write_status(status_file,"LOAD_PACKAGES", "loading")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(cmdstanr)
  library(posterior)
  library(tibble)
  library(shigella)
})

pkg_versions <- c(
  R            = R.version.string,
  cmdstanr     = as.character(packageVersion("cmdstanr")),
  posterior    = as.character(packageVersion("posterior")),
  shigella     = as.character(packageVersion("shigella")),
  serodynamics = as.character(packageVersion("serodynamics"))
)
cmdstan_ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) "UNKNOWN")
cat("=== Package versions ===\n")
for (n in names(pkg_versions)) cat(sprintf("  %-15s %s\n", n, pkg_versions[n]))
cat(sprintf("  %-15s %s\n\n", "cmdstan", cmdstan_ver))
write_status(status_file,"LOAD_PACKAGES", "OK")

# ----- 2. Compile dir — SLURM-specific (per-task subdir) -----
write_status(status_file,"COMPILE_DIR", "setting up")
compile_dir <- Sys.getenv("STAN_COMPILE_DIR", unset = "")
if (compile_dir == "") {
  # Fallback if sbatch didn't set it
  user <- Sys.getenv("USER", unset = "unknown")
  compile_dir <- file.path("/tmp", user, sprintf("cmdstan_bin_phase1_%s", job_id))
}
if (!dir.exists(compile_dir)) {
  dir.create(compile_dir, recursive = TRUE, mode = "0755")
}
cat(sprintf("=== compile_dir: %s ===\n", compile_dir))
cat(sprintf("    existing files: %d\n\n", length(list.files(compile_dir))))
write_status(status_file,"COMPILE_DIR", "OK")

# ----- 3. Simulate (IDENTICAL to Phase 0 — same seed, same n) -----
write_status(status_file,"SIMULATE", "running")
set.seed(20260513)
true_rho_B   <- 0.6
omega_B_true <- matrix(c(1, true_rho_B, true_rho_B, 1), 2, 2)

sim_data <- sim_correlated_case_data(
  n                 = 5,
  omega_B           = omega_B_true,
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)
cat(sprintf("=== Simulated: n=%d, rho_B=%.1f, total rows=%d ===\n\n",
            length(unique(sim_data$id)), true_rho_B, nrow(sim_data)))
write_status(status_file,"SIMULATE", "OK")

# ----- 4. Save STARTED placeholder so crash mid-fit is still informative -----
out_file <- sprintf("outputs/phase1/one_fit_n5_jobid_%s.rds", job_id)
saveRDS(
  list(
    scenario   = "phase1_slurm_single",
    status     = "FIT_STARTED",
    job_id     = job_id,
    started_at = format(Sys.time()),
    slurm_env  = slurm_env,
    pkg_versions = pkg_versions,
    cmdstan_version = cmdstan_ver,
    true_rho_B = true_rho_B,
    n = 5
  ),
  out_file
)

# ----- 5. Run fit -----
write_status(status_file,"FIT", "running")
cat("=== Fitting model_2 (Kronecker) ===\n")
t_start <- Sys.time()

fit <- tryCatch({
  run_mod_stan(
    data            = sim_data,
    model           = "model_2",
    chains          = 2,
    iter_warmup     = 500,
    iter_sampling   = 500,
    parallel_chains = 2,
    adapt_delta     = 0.95,
    max_treedepth   = 12,
    init            = 0.1,
    with_post       = TRUE,
    stan_dir        = "inst/stan",
    compile_dir     = compile_dir,
    refresh         = 100,
    show_messages   = TRUE
  )
}, error = function(e) {
  cat("\n  [FIT ERROR]:", conditionMessage(e), "\n")
  cat("  [STACK TRACE]:\n")
  print(sys.calls())
  write_status(status_file,"FIT", paste("CRASHED:", conditionMessage(e)))
  saveRDS(
    list(scenario = "phase1_slurm_single", status = "FIT_FAILED",
         job_id = job_id,
         error = conditionMessage(e),
         crashed_at = format(Sys.time()),
         slurm_env = slurm_env,
         pkg_versions = pkg_versions),
    out_file
  )
  NULL
})

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\n=== Fit elapsed: %.2f min ===\n\n", elapsed))

if (is.null(fit)) {
  cat(strrep("!", 70), "\n", sep = "")
  cat(" PHASE 1 RESULT: FIT CRASHED INSIDE SLURM\n")
  cat(" Compare with outputs/phase0/one_fit_n5.rds to determine:\n")
  cat("   - If Phase 0 OK but Phase 1 FAIL → Slurm env issue\n")
  cat("   - If both fail → code / model identifiability issue\n")
  cat(strrep("!", 70), "\n", sep = "")
  write_status(status_file,"DONE", "Phase 1 FAILED")
  quit(status = 1)
}

write_status(status_file,"FIT", "OK")

# ----- 6. Diagnostics -----
write_status(status_file,"DIAG", "extracting")
sf <- attr(fit, "stan_fit")[[1]]
diag <- sf$diagnostic_summary(diagnostics = c("divergences", "treedepth", "ebfmi"))
total_iters <- 2 * 500

cat("=== Diagnostics ===\n")
cat(sprintf("  divergent:        %d / %d (%.2f%%)\n",
            sum(diag$num_divergent), total_iters,
            100 * sum(diag$num_divergent) / total_iters))
cat(sprintf("  max-treedepth:    %d / %d (%.2f%%)\n",
            sum(diag$num_max_treedepth), total_iters,
            100 * sum(diag$num_max_treedepth) / total_iters))
cat(sprintf("  E-BFMI:           %s\n",
            paste(sprintf("%.3f", diag$ebfmi), collapse = ", ")))

draws_summary <- posterior::summarise_draws(
  sf$draws(variables = "Omega_B[1,2]"),
  "median", "mean", "sd",
  ~quantile(.x, c(0.025, 0.975), na.rm = TRUE),
  "ess_bulk", "rhat"
)
cat("\n  Omega_B[1,2] summary:\n")
print(draws_summary)
write_status(status_file,"DIAG", "OK")

# ----- 7. Save bundle -----
write_status(status_file,"SAVE", "writing rds")
result_bundle <- list(
  scenario            = "phase1_slurm_single",
  status              = "OK",
  job_id              = job_id,
  elapsed_min         = elapsed,
  started_at          = format(t_start),
  completed_at        = format(Sys.time()),
  slurm_env           = slurm_env,
  pkg_versions        = pkg_versions,
  cmdstan_version     = cmdstan_ver,
  true_rho_B          = true_rho_B,
  n_subjects          = 5,
  fit_settings        = list(chains = 2, warmup = 500, sampling = 500,
                             adapt_delta = 0.95, max_treedepth = 12),
  diagnostic_summary  = diag,
  omega_B_summary     = draws_summary,
  rho_B_posterior     = as.vector(posterior::as_draws_array(
                          sf$draws("Omega_B[1,2]")))
)
saveRDS(result_bundle, out_file)

# ----- 8. Compare with Phase 0 result if it exists -----
phase0_file <- "outputs/phase0/one_fit_n5.rds"
if (file.exists(phase0_file)) {
  ph0 <- readRDS(phase0_file)
  if (!is.null(ph0$omega_B_summary)) {
    cat("\n=== Phase 0 vs Phase 1 comparison ===\n")
    cmp <- data.frame(
      metric            = c("status", "elapsed_min", "post_median",
                            "post_lo_2.5", "post_hi_97.5",
                            "ess_bulk", "rhat",
                            "n_divergent", "n_treedepth"),
      phase0 = c(ph0$status,
                 round(ph0$elapsed_min, 2),
                 round(ph0$omega_B_summary$median, 3),
                 round(ph0$omega_B_summary$`2.5%`, 3),
                 round(ph0$omega_B_summary$`97.5%`, 3),
                 round(ph0$omega_B_summary$ess_bulk, 0),
                 round(ph0$omega_B_summary$rhat, 3),
                 sum(ph0$diagnostic_summary$num_divergent),
                 sum(ph0$diagnostic_summary$num_max_treedepth)),
      phase1 = c("OK",
                 round(elapsed, 2),
                 round(draws_summary$median, 3),
                 round(draws_summary$`2.5%`, 3),
                 round(draws_summary$`97.5%`, 3),
                 round(draws_summary$ess_bulk, 0),
                 round(draws_summary$rhat, 3),
                 sum(diag$num_divergent),
                 sum(diag$num_max_treedepth))
    )
    print(cmp, row.names = FALSE)

    
    saveRDS(cmp, sprintf("outputs/phase1/p0_vs_p1_comparison_%s.rds", job_id))
  }
} else {
  cat("\n  [INFO] Phase 0 result not found — comparison skipped.\n")
  cat("         Run phase0_interactive_reproducibility.R first if you want\n")
  cat("         direct apples-to-apples comparison.\n")
}

write_status(status_file,"SAVE", "OK")
cat("\n=== Phase 1 complete ===\n")
cat(sprintf("    Results: %s\n\n", out_file))
write_status(status_file,"DONE", "Phase 1 OK")
