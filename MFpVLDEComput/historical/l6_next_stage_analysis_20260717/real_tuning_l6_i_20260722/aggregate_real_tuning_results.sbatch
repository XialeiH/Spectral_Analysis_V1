#!/bin/bash
#SBATCH --job-name=real_tune_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${REAL_TUNING_CODE_ROOT}','-begin'); aggregate_real_tuning_results"
