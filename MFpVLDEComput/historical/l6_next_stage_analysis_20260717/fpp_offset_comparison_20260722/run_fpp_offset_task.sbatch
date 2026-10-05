#!/bin/bash
#SBATCH --job-name=fpp_offset
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${FPP_OFFSET_MAIN_ROOT}'); addpath('${FPP_OFFSET_MAIN_ROOT}/Utils'); addpath('${FPP_OFFSET_H96_RUNTIME}','-begin'); addpath('${FPP_OFFSET_CODE_ROOT}','-begin'); run_fpp_offset_task"
