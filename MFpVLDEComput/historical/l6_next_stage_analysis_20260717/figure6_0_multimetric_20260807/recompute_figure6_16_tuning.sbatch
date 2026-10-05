#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-8%9
#SBATCH --job-name=f616fix

set -euo pipefail
run_root="${FIG60_RUN_ROOT:?FIG60_RUN_ROOT is required}"
code_root="${run_root}/code"
export FIG616_DATASET_FILE="${run_root}/results/contrast_100/figure6_0_contrast100_dataset.mat"
export FIG616_CORRECTED_ROOT="${run_root}/results/figure6_16_tuning_corrected"
export FIG616_SELECTION_INDEX="${SLURM_ARRAY_TASK_ID}"

module purge
module load matlab/2025b
matlab -batch "maxNumCompThreads(4); addpath('${FIG60_MAIN_ROOT}'); addpath('${FIG60_MAIN_ROOT}/Utils'); addpath('${FIG60_H96_RUNTIME}','-begin'); addpath('${code_root}','-begin'); recompute_figure6_16_tuning"
