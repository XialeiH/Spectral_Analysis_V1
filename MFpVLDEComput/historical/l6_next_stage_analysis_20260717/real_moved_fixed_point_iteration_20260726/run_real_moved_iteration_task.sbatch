#!/bin/bash
#SBATCH --job-name=real_move_iter
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=03:30:00

set -euo pipefail
module purge
module load matlab/2025b

mkdir -p "${REAL_MOVED_OUTPUT_ROOT}/logs"
matlab -batch "addpath('${REAL_MOVED_MAIN_ROOT}'); addpath('${REAL_MOVED_MAIN_ROOT}/Utils'); addpath('${REAL_MOVED_H96_RUNTIME}','-begin'); addpath('${REAL_MOVED_CODE_ROOT}','-begin'); run_real_moved_iteration_task"
