#!/bin/bash
#SBATCH --job-name=tshc_gap
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail
module purge
module load matlab/2025b
mkdir -p "${TSHC_OUTPUT_ROOT}"
matlab -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); addpath('${TSHC_CODE_ROOT}','-begin'); check_baseline_singular_gap"
