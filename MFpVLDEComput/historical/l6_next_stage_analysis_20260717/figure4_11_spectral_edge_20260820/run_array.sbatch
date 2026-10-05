#!/bin/bash
#SBATCH --job-name=f4_11_edge
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --array=1-16%16
#SBATCH --output=logs/edge_%A_%a.out
#SBATCH --error=logs/edge_%A_%a.err

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${RUN_CODE}'); run_spectral_edge_condition(${SLURM_ARRAY_TASK_ID},'${GEOMETRY_ROOT}','${RESULT_ROOT}','${CODE_ROOT}');"
