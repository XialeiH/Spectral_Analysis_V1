#!/bin/bash
#SBATCH --job-name=l6i_selected_fp
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${FULL_GRID_MAIN_ROOT}'); addpath('${FULL_GRID_MAIN_ROOT}/Utils'); addpath('${FULL_GRID_H96_RUNTIME}','-begin'); addpath('${FULL_GRID_CODE_ROOT}','-begin'); addpath('${SELECTED_CODE_ROOT}','-begin'); plot_selected_fixed_points('${FULL_GRID_SETUP_FILE}','${SELECTED_OUTPUT_ROOT}')"
