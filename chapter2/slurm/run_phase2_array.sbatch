#!/bin/bash
# ==========================================================================
# run_phase2_array.sbatch — SLURM array job for Chapter 2 Phase 2
#
# Submits R*3 array tasks (R=500 -> 1500 tasks) to run Phase 2 simulation
# in parallel. Each task runs ONE replicate of one scenario.
#
# Submit with:
#   cd ~/chapter2
#   sbatch slurm/run_phase2_array.sbatch
#
# Monitor with:
#   squeue -u $USER
#   sacct -j <job_id>
#   tail -f logs/phase2_<job_id>_<task_id>.out
#
# Aggregate results after all tasks complete:
#   Rscript scripts/02_collect_array.R
# ==========================================================================

#SBATCH --job-name=ch2_phase2
#SBATCH --output=logs/phase2_%A_%a.out
#SBATCH --error=logs/phase2_%A_%a.err

# Array configuration:
#   R = 500 reps per scenario
#   3 scenarios (A, B, C)
#   = 1500 array tasks total
# Adjust --array as needed, e.g. --array=1-1500%50 to limit concurrency
#SBATCH --array=1-1500%50

# Resource per task:
#   4 chains x ~3000 iterations -> ~30-60 min per task on Shiva
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# ----- Optional: partition / qos (uncomment if Shiva requires) -----
# #SBATCH --partition=high
# #SBATCH --qos=high

set -euo pipefail

# ==========================================================================
# 1. Conda environment activation
# ==========================================================================
# Initialize conda for non-interactive shells
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate r_chapter2

echo "=== Job: $SLURM_JOB_NAME  Task: $SLURM_ARRAY_TASK_ID ==="
echo "=== Host: $(hostname)  Date: $(date) ==="
echo "=== R: $(which R)  Rscript: $(which Rscript) ==="

# ==========================================================================
# 2. Configure compile output dir (writable + executable on Shiva)
# ==========================================================================
# Per-user /tmp dir to avoid collisions between concurrent tasks
export STAN_COMPILE_DIR="/tmp/$USER/cmdstan_bin"
mkdir -p "$STAN_COMPILE_DIR"

echo "=== Stan compile dir: $STAN_COMPILE_DIR ==="

# ==========================================================================
# 3. Configure thread count for nested parallelism inside R
# ==========================================================================
# We use 4 chains in parallel; each chain on 1 core. Disable BLAS threading
# to avoid oversubscription.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ==========================================================================
# 4. R_TOTAL — must match scripts/02_run_array.R expectation
# ==========================================================================
export R_TOTAL=500

# ==========================================================================
# 5. Run the R script for this task
# ==========================================================================
cd "$HOME/chapter2"

mkdir -p logs outputs/02_array

Rscript scripts/02_run_array.R

echo "=== Task $SLURM_ARRAY_TASK_ID finished at $(date) ==="
