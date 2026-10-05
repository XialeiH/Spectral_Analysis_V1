#!/bin/bash
#SBATCH --job-name=l6_ei_hc_grid
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
mkdir -p "${L6EI_OUTPUT_ROOT}/logs"
task_id="${SLURM_ARRAY_TASK_ID:-${L6EI_TASK_ID:-unknown}}"
temporary_log="${L6EI_OUTPUT_ROOT}/logs/task_${SLURM_ARRAY_JOB_ID:-single}_${task_id}.tmp"
error_log="${temporary_log%.tmp}.err"
if matlab -batch "addpath('${L6EI_MAIN_ROOT}'); addpath('${L6EI_MAIN_ROOT}/Utils'); addpath('${L6EI_H96_RUNTIME}','-begin'); addpath('${L6EI_CODE_ROOT}','-begin'); run_l6_ei_grid_task" >"$temporary_log" 2>&1; then
  rm -f "$temporary_log"
else
  mv "$temporary_log" "$error_log"
  exit 1
fi
