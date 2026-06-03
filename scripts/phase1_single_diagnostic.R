# ==========================================================================
# phase1_single_diagnostic.R  —  n = 5 pilot
#
# Launched by slurm/phase1_single.sbatch
#
# Purpose: fit model_2 inside a SLURM sbatch job and compare with the
# Phase 0 interactive baseline to isolate SLURM-vs-code attribution.
# ==========================================================================
suppressPackageStartupMessages(library(shigella))

shigella:::run_phase1_diagnostic(
  n             = 5,
  iter_warmup   = 500,
  iter_sampling = 500,
  tag           = "n5"
)
