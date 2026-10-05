#!/bin/bash
#SBATCH --job-name=global_offset
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${GLOBAL_OFFSET_MAIN_ROOT}'); addpath('${GLOBAL_OFFSET_MAIN_ROOT}/Utils'); addpath('${GLOBAL_OFFSET_H96_RUNTIME}','-begin'); addpath('${GLOBAL_OFFSET_CODE_ROOT}','-begin'); run_global_offset_task"
