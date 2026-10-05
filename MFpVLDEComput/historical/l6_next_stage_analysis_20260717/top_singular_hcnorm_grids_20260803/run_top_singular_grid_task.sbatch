#!/bin/bash
#SBATCH --job-name=tshc_grid
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:30:00

set -euo pipefail
module purge
module load matlab/2025b
temporary_log="${TSHC_OUTPUT_ROOT}/logs/task_${SLURM_ARRAY_JOB_ID:-single}_${SLURM_ARRAY_TASK_ID:-${TSHC_TASK_ID:-unknown}}.tmp"
error_log="${temporary_log%.tmp}.err"
if matlab -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); addpath('${TSHC_MAIN_ROOT}'); addpath('${TSHC_MAIN_ROOT}/Utils'); addpath('${TSHC_H96_RUNTIME}','-begin'); addpath('${TSHC_CODE_ROOT}','-begin'); run_top_singular_grid_task" >"$temporary_log" 2>&1; then
  rm -f "$temporary_log"
else
  mv "$temporary_log" "$error_log"
  exit 1
fi
