#!/bin/bash
#SBATCH --job-name=l6_i_offset
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${OFFSET_MAIN_ROOT}'); addpath('${OFFSET_MAIN_ROOT}/Utils'); addpath('${OFFSET_H96_RUNTIME}','-begin'); addpath('${OFFSET_CODE_ROOT}','-begin'); run_l6_i_offset_task"
