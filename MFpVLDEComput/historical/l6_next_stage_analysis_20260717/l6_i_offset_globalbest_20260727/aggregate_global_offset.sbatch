#!/bin/bash
#SBATCH --job-name=global_offset_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${GLOBAL_OFFSET_CODE_ROOT}','-begin'); aggregate_global_offset"
