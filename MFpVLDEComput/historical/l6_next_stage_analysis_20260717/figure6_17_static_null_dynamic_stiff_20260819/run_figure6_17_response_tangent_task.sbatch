#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=1-28%28
#SBATCH --job-name=f617T

set -euo pipefail
module purge
module load matlab/2025b
mkdir -p "${FIG617_RUN_ROOT}/results/figure6_17/logs"
matlab -batch "maxNumCompThreads(4); addpath('${FIG60_MAIN_ROOT}'); addpath('${FIG60_MAIN_ROOT}/Utils'); addpath('${FIG60_H96_RUNTIME}','-begin'); addpath('${FIG60_CODE_ROOT}','-begin'); run_figure6_17_response_tangent_task" \
  >"${FIG617_RUN_ROOT}/results/figure6_17/logs/task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
  2>"${FIG617_RUN_ROOT}/results/figure6_17/logs/task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
