#!/usr/bin/env bash
#SBATCH --job-name=ji_eigenspace_angles
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=1-4%4

set -euo pipefail

: "${JI_CODE_ROOT:?JI_CODE_ROOT is required}"
: "${JI_SOURCE_ROOT:?JI_SOURCE_ROOT is required}"
: "${JI_OUTPUT_ROOT:?JI_OUTPUT_ROOT is required}"

module purge
module load matlab/2025b

matlab -batch "addpath('${JI_CODE_ROOT}'); run_ji_top_real_eigenspace_angles_task([], '${JI_SOURCE_ROOT}', '${JI_OUTPUT_ROOT}');"
