#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-100%25
#SBATCH --job-name=fig1g_timing

set -euo pipefail
RUN_ROOT="$1"
module purge
module load matlab/2025b
mkdir -p "$RUN_ROOT/results" "$RUN_ROOT/logs"
matlab -batch "addpath('$RUN_ROOT/code'); run_figure1g_trial('$RUN_ROOT', '$SLURM_ARRAY_TASK_ID');"
