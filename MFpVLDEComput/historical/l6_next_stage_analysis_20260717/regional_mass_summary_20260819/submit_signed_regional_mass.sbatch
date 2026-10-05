#!/bin/bash
#SBATCH --job-name=regional_mass
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-16%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -euo pipefail
module purge
module load matlab/2025b

RUN_ROOT=${1:?run root is required}
SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/l6_mechanism_derivative_clamp_20260628_033413
mkdir -p "${RUN_ROOT}/results" "${RUN_ROOT}/logs"
matlab -batch "addpath('${RUN_ROOT}/code'); compute_signed_regional_mass_task(${SLURM_ARRAY_TASK_ID}, '${SOURCE_ROOT}', '${RUN_ROOT}/results');"
