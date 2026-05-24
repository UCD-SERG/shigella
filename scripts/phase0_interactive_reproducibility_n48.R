# ============================================================================
# phase0_interactive_reproducibility_n48.R  —  n = 48 full cohort
#
# Execution:
#   salloc --time=08:00:00 --cpus-per-task=2 --mem=20G
#   Rscript scripts/phase0_interactive_reproducibility_n48.R
#   exit
#
# Purpose: reproduce a Phase 1 (sbatch) fit under interactive SLURM to confirm
# determinism across allocation modes. n=48 version (full cohort size).
# Run phase0_interactive_reproducibility.R (n=5) first to confirm the
# pipeline works before committing to this longer run.
# ============================================================================
setwd("~/shigella")
suppressPackageStartupMessages(library(shigella))

shigella:::run_phase0_diagnostic(
  n             = 48,
  iter_warmup   = 1000,
  iter_sampling = 1000,
  tag           = "n48"
)
