#!/bin/bash
#SBATCH --job-name=f4D_plot
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${RUN_CODE}'); plot_finite_time_gain('${RESULT_ROOT}','${FIGURE_ROOT}');"
