#!/bin/bash
#SBATCH --job-name=real_tune
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
matlab -batch "addpath('${REAL_TUNING_MAIN_ROOT}'); addpath('${REAL_TUNING_MAIN_ROOT}/Utils'); addpath('${REAL_TUNING_H96_RUNTIME}','-begin'); addpath('${REAL_TUNING_CODE_ROOT}','-begin'); run_real_tuning_task"
