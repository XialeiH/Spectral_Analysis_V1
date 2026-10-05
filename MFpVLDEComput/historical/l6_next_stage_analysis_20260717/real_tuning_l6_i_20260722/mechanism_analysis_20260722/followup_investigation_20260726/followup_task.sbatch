#!/bin/bash
#SBATCH --job-name=h96_followup
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=18G
#SBATCH --time=03:45:00

set -euo pipefail
module purge
module load matlab/2025b
mkdir -p "${FOLLOWUP_OUTPUT_ROOT}/logs"
matlab -batch "addpath('${FOLLOWUP_MAIN_ROOT}'); addpath('${FOLLOWUP_MAIN_ROOT}/Utils'); addpath('${FOLLOWUP_H96_RUNTIME}','-begin'); addpath('${FOLLOWUP_CODE_ROOT}','-begin'); ${FOLLOWUP_ENTRY}"
