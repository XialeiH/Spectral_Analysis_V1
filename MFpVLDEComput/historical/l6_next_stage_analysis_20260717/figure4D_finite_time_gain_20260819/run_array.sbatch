#!/bin/bash
#SBATCH --job-name=f4D_gain
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --array=1-16%16

set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('${RUN_CODE}'); run_gain_condition(${SLURM_ARRAY_TASK_ID},'${GEOMETRY_ROOT}','${RESULT_ROOT}','${SHARED_CODE}');"
