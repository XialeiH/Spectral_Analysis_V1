#!/bin/bash
#SBATCH --job-name=real_iter_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail
module purge
module load matlab/2025b

matlab -batch "addpath('${REAL_MOVED_MAIN_ROOT}'); addpath('${REAL_MOVED_MAIN_ROOT}/Utils'); addpath('${REAL_MOVED_H96_RUNTIME}','-begin'); addpath('${REAL_MOVED_CODE_ROOT}','-begin'); aggregate_real_moved_iteration_results"
