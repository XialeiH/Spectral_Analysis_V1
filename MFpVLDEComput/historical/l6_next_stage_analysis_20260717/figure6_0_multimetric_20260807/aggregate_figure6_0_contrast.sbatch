#!/bin/bash
#SBATCH --job-name=fig60_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${FIG60_CODE_ROOT}','-begin'); aggregate_figure6_0_contrast"
