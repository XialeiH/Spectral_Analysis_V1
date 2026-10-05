#!/bin/bash
#SBATCH --job-name=pair_hc_ode
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b
mkdir -p "${PAIR_GRID_ODE_OUTPUT_ROOT}/logs"
temporary_log="${PAIR_GRID_ODE_OUTPUT_ROOT}/logs/task_${SLURM_ARRAY_JOB_ID:-single}_${SLURM_ARRAY_TASK_ID:-${PAIR_GRID_ODE_INDEX:-unknown}}.tmp"
error_log="${temporary_log%.tmp}.err"
if matlab -batch "addpath('${PAIR_GRID_MAIN_ROOT}'); addpath('${PAIR_GRID_MAIN_ROOT}/Utils'); addpath('${PAIR_GRID_H96_RUNTIME}','-begin'); addpath('${PAIR_GRID_CODE_ROOT}','-begin'); run_pair_grid_ode_fallback_task" >"$temporary_log" 2>&1; then
  rm -f "$temporary_log"
else
  mv "$temporary_log" "$error_log"
  exit 1
fi
