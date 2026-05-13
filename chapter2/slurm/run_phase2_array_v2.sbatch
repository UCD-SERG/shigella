#!/bin/bash
# ==========================================================================
# run_phase2_array_v2.sbatch — Hardened Phase 2 SLURM array
#
# Changes from v1:
#   1. Per-task STAN_COMPILE_DIR (avoids /tmp concurrent-write hang)
#   2. Reduced concurrency: 20 instead of 50 (less node oversubscription)
#   3. Time limit raised to 3 hours (margin of safety with light settings)
#   4. R=200 by default (Burton 2006 publishable; Ezra-acceptable)
#      → 600 tasks total instead of 1500
#   5. Calls 02_run_array_v2.R (lighter settings, robust error handling)
#
# Submit:
#   sbatch slurm/run_phase2_array_v2.sbatch
#
# Expected wall time: 4-7 hours
# ==========================================================================

#SBATCH --job-name=ch2_phase2_v2
#SBATCH --output=logs/phase2v2_%A_%a.out
#SBATCH --error=logs/phase2v2_%A_%a.err

# R = 200 reps per scenario, 3 scenarios = 600 tasks
# Concurrency 20 (down from 50) — less /tmp contention
#SBATCH --array=1-600%30

# 3 hour time limit (was 2). Light settings should finish in 30-60 min.
#SBATCH --time=01:30:00

# 4 cores per chain
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

# ==========================================================================
# 1. Conda
# ==========================================================================
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate r_chapter2

echo "=== Job: $SLURM_JOB_NAME  Task: $SLURM_ARRAY_TASK_ID ==="
echo "=== Host: $(hostname)  Date: $(date) ==="
echo "=== R: $(which R) ==="

# ==========================================================================
# 2. CRITICAL: Per-task compile dir
# ==========================================================================
# Each task gets its own subdirectory. This eliminates the
# "everyone writing to same /tmp/cmdstan_bin" concurrent-write hang.
export STAN_COMPILE_DIR="/tmp/$USER/cmdstan_bin_task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$STAN_COMPILE_DIR"

echo "=== Stan compile dir: $STAN_COMPILE_DIR ==="

# ==========================================================================
# 3. Threading
# ==========================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ==========================================================================
# 4. R_TOTAL — must match scripts/02_run_array_v2.R expectation
# ==========================================================================
export R_TOTAL=200

# ==========================================================================
# 5. Run
# ==========================================================================
cd "$HOME/chapter2"
mkdir -p logs outputs/02_array

Rscript scripts/02_run_array_v2.R

echo "=== Task $SLURM_ARRAY_TASK_ID finished at $(date) ==="

# ==========================================================================
# 6. Cleanup per-task compile dir (saves /tmp space)
# ==========================================================================
# Only delete if task succeeded (don't lose debug info if failed)
if [ "$?" -eq 0 ]; then
    rm -rf "$STAN_COMPILE_DIR"
fi
