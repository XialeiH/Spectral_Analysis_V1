#!/bin/bash
#SBATCH --job-name=l6i_hc_grid
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:45:00

set -euo pipefail
module purge
module load matlab/2025b
temporary_log="${FULL_GRID_OUTPUT_ROOT}/logs/task_${SLURM_ARRAY_JOB_ID:-single}_${SLURM_ARRAY_TASK_ID:-${FULL_GRID_TASK_ID:-unknown}}.tmp"
error_log="${temporary_log%.tmp}.err"
if matlab -batch "addpath('${FULL_GRID_MAIN_ROOT}'); addpath('${FULL_GRID_MAIN_ROOT}/Utils'); addpath('${FULL_GRID_H96_RUNTIME}','-begin'); addpath('${FULL_GRID_CODE_ROOT}','-begin'); run_full_grid_task" >"$temporary_log" 2>&1; then
  rm -f "$temporary_log"
else
  mv "$temporary_log" "$error_log"
  exit 1
fi
