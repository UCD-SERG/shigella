# ==========================================================================
# 02_run_array_v2.R — Hardened SLURM array task version
#
# Changes from v1:
#   1. LIGHTER settings (500 warmup + 500 sampling, adapt_delta = 0.92)
#   2. Per-task STAN_COMPILE_DIR support (avoids /tmp concurrent-write hang)
#   3. OUTPUT_DIR env var support (for test runs vs full runs)
#   4. Much more verbose logging at every step
#   5. Saves a status file IMMEDIATELY at start so we know task started
#   6. Robust error handling — write FAILED status even if R crashes
# ==========================================================================

setwd("~/chapter2")

cat("\n=== 02_run_array_v2.R START ===\n")
cat("Time:", format(Sys.time()), "\n")

# Track what step we're at — so even if we crash mid-way, the log shows where
.STEP <- function(msg) cat(sprintf("[STEP] %s\n", msg))

.STEP("Load packages")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(serodynamics)
  library(cmdstanr)
  library(posterior)
  library(cli)
  library(tibble)
})

.STEP("Source helpers")
source("R/prep_data_stan.R")
source("R/prep_priors_stan.R")
source("R/postprocess_stan_output.R")
source("R/run_mod_stan.R")
source("R/sim_correlated_case_data.R")

# ==========================================================================
# 1. Read SLURM array task ID + R_TOTAL
# ==========================================================================
.STEP("Read SLURM env")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
R_TOTAL <- as.integer(Sys.getenv("R_TOTAL", unset = "200"))

cat(sprintf("  task_id = %d\n", task_id))
cat(sprintf("  R_TOTAL = %d\n", R_TOTAL))

if (is.na(task_id) || task_id < 1) {
  stop("Invalid SLURM_ARRAY_TASK_ID: ", task_id)
}

# ==========================================================================
# 2. Decode task ID -> (scenario, rep)
# ==========================================================================
.STEP("Decode task ID")
scenario_idx <- ((task_id - 1) %/% R_TOTAL) + 1L  # 1, 2, or 3
rep_idx      <- ((task_id - 1) %%  R_TOTAL) + 1L  # 1..R

scenarios <- list(
  list(name = "A", n = 48, rho_B = 0.6),
  list(name = "B", n = 11, rho_B = 0.6),
  list(name = "C", n = 48, rho_B = 0.0)
)

if (scenario_idx > length(scenarios)) {
  stop("scenario_idx out of range: ", scenario_idx)
}
scn <- scenarios[[scenario_idx]]

cat(sprintf("  Scenario %s, rep %d (n=%d, true rho_B=%.1f)\n",
            scn$name, rep_idx, scn$n, scn$rho_B))

# ==========================================================================
# 3. Output path
# ==========================================================================
.STEP("Set output path")
out_dir <- Sys.getenv("OUTPUT_DIR", unset = "outputs/02_array")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("scenario_%s_rep_%04d.rds", scn$name, rep_idx))

cat(sprintf("  Output: %s\n", out_file))

# Skip if already done
if (file.exists(out_file)) {
  cat("  [SKIP] Already done\n")
  quit(status = 0)
}

# Write a placeholder so we know the task started even if it dies later
saveRDS(list(scenario = scn$name, rep = rep_idx, task_id = task_id,
             status = "STARTED", started_at = format(Sys.time())),
        out_file)

# ==========================================================================
# 4. Compile dir — CRITICAL for SLURM array (per-task subdir)
# ==========================================================================
.STEP("Set up compile dir")
compile_dir <- Sys.getenv("STAN_COMPILE_DIR", unset = "")
if (compile_dir == "") {
  user <- Sys.getenv("USER", unset = "default")
  compile_dir <- file.path("/tmp", user, "cmdstan_bin")
}
if (!dir.exists(compile_dir)) {
  dir.create(compile_dir, recursive = TRUE, mode = "0755")
}
cat(sprintf("  compile_dir = %s\n", compile_dir))
cat(sprintf("  Existing files: %d\n", length(list.files(compile_dir))))

# ==========================================================================
# 5. Run one fit
# ==========================================================================
make_omega_2x2 <- function(rho) {
  matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
}

.STEP("Begin fit")
t0 <- Sys.time()

set.seed(2026 * 1000 + scenario_idx * R_TOTAL + rep_idx)

Omega_B_true <- make_omega_2x2(scn$rho_B)

# ----- 5a. Simulate -----
.STEP("Simulate data")
sim_dat <- tryCatch({
  sim_correlated_case_data(
    n                 = scn$n,
    Omega_B           = Omega_B_true,
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )
}, error = function(e) {
  cat("  SIM ERROR:", conditionMessage(e), "\n")
  NULL
})

if (is.null(sim_dat)) {
  result <- list(scenario = scn$name, rep = rep_idx,
                 task_id = task_id, status = "SIM_FAILED")
  saveRDS(result, out_file)
  quit(status = 0)
}
cat(sprintf("  sim_dat rows: %d\n", nrow(sim_dat)))

# ----- 5b. Fit -----
.STEP("Fit (run_mod_stan)")
fit <- tryCatch({
  run_mod_stan(
    data            = sim_dat,
    model           = "model_2",
    chains          = 4,
    iter_warmup     = 300,    # was 1500 — much lighter
    iter_sampling   = 400,    # was 1500
    parallel_chains = 4,
    adapt_delta     = 0.90,   # was 0.99 — faster sampling
    max_treedepth   = 11,     # was 12
    init            = 0.1,
    with_post       = TRUE,
    stan_dir        = "inst/stan",
    compile_dir     = compile_dir,
    refresh         = 100,
    show_messages   = FALSE
  )
}, error = function(e) {
  cat("  FIT ERROR:", conditionMessage(e), "\n")
  NULL
})

if (is.null(fit)) {
  result <- list(scenario = scn$name, rep = rep_idx,
                 task_id = task_id, status = "FIT_FAILED",
                 elapsed_min = as.numeric(Sys.time() - t0, units = "mins"))
  saveRDS(result, out_file)
  quit(status = 0)
}

# ----- 5c. Extract posterior -----
.STEP("Extract Omega_B posterior")
rho_B_post <- tryCatch({
  sf <- attr(fit, "stan_fit")[[1]]
  omega_B_draws <- posterior::as_draws_df(
    sf$draws(variables = "Omega_B[1,2]")
  )
  omega_B_draws[["Omega_B[1,2]"]]
}, error = function(e) {
  cat("  EXTRACT ERROR:", conditionMessage(e), "\n")
  NULL
})

if (is.null(rho_B_post)) {
  result <- list(scenario = scn$name, rep = rep_idx,
                 task_id = task_id, status = "EXTRACT_FAILED",
                 elapsed_min = as.numeric(Sys.time() - t0, units = "mins"))
  saveRDS(result, out_file)
  quit(status = 0)
}

# ----- 5d. Diagnostics + save -----
.STEP("Save result")
n_divergent <- tryCatch({
  sf <- attr(fit, "stan_fit")[[1]]
  sum(sf$diagnostic_summary()$num_divergent)
}, error = function(e) NA_integer_)

elapsed <- as.numeric(Sys.time() - t0, units = "mins")

result <- list(
  scenario          = scn$name,
  rep               = rep_idx,
  task_id           = task_id,
  status            = "OK",
  true_rho_B        = scn$rho_B,
  est_rho_B_median  = median(rho_B_post),
  est_rho_B_mean    = mean(rho_B_post),
  est_rho_B_lo      = quantile(rho_B_post, 0.025, names = FALSE),
  est_rho_B_hi      = quantile(rho_B_post, 0.975, names = FALSE),
  bias              = median(rho_B_post) - scn$rho_B,
  n_divergent       = n_divergent,
  elapsed_min       = elapsed,
  n_post_draws      = length(rho_B_post)
)

saveRDS(result, out_file)

cat(sprintf("\n  rho_B = %.3f [%.3f, %.3f], bias = %+.3f, %d div, %.1f min\n",
            result$est_rho_B_median, result$est_rho_B_lo, result$est_rho_B_hi,
            result$bias, n_divergent, elapsed))

cat(sprintf("\n=== Task %d DONE: %s ===\n", task_id, out_file))
