#!/bin/bash
#SBATCH --job-name=pair_hc_smoke
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${PAIR_GRID_MAIN_ROOT}'); addpath('${PAIR_GRID_MAIN_ROOT}/Utils'); addpath('${PAIR_GRID_H96_RUNTIME}','-begin'); addpath('${PAIR_GRID_CODE_ROOT}','-begin'); run_pair_grid_smoke"
