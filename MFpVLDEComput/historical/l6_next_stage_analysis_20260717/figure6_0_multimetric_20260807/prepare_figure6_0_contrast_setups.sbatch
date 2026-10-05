#!/bin/bash
#SBATCH --job-name=fig60_prep
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); addpath('${FIG60_MAIN_ROOT}'); addpath('${FIG60_MAIN_ROOT}/Utils'); addpath('${FIG60_H96_RUNTIME}','-begin'); addpath('${FIG60_CODE_ROOT}','-begin'); prepare_figure6_0_contrast_setups"
