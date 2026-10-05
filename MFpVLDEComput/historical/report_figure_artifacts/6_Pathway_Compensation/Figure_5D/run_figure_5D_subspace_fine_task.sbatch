#!/bin/bash
#SBATCH --job-name=fig5d_subspace_fine
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:30:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${FIGURE5D_MAIN_ROOT}'); addpath('${FIGURE5D_MAIN_ROOT}/Utils'); addpath('${FIGURE5D_H96_RUNTIME}','-begin'); addpath('${FIGURE5D_CODE_ROOT}','-begin'); run_figure_5D_subspace_fine_task"
