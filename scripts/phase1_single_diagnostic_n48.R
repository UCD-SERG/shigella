# ==========================================================================
# phase1_single_diagnostic_n48.R  —  n = 48 full cohort
#
# Launched by slurm/phase1_single_n48.sbatch
#
# Purpose: fit model_2 inside a SLURM sbatch job and compare with the
# Phase 0 interactive baseline to isolate SLURM-vs-code attribution.
# n=48 version (full cohort size).
# ==========================================================================
suppressPackageStartupMessages(library(shigella))

run_phase1_diagnostic(
  n             = 48,
  iter_warmup   = 1000,
  iter_sampling = 1000,
  tag           = "n48"
)
