#!/bin/bash
#SBATCH --job-name=fpp_ode
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:30:00

set -euo pipefail
module purge
module load matlab/2025b

mkdir -p "${FPP_OUTPUT_ROOT}/logs"
matlab -batch "addpath('${FPP_MAIN_ROOT}'); addpath('${FPP_MAIN_ROOT}/Utils'); addpath('${H96_RUNTIME_DIR}','-begin'); addpath('${FPP_CODE_ROOT}','-begin'); addpath('${FPP_CODE_ROOT}/Utils','-begin'); run_fpp_ode_task"
