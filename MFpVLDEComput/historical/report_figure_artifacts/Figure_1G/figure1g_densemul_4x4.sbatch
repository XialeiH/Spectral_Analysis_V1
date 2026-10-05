#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --array=1-100%6
#SBATCH --job-name=fig1g_d4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$RUN_ROOT/code'); run_figure1g_field_trial_densemul('$RUN_ROOT',4,$SLURM_ARRAY_TASK_ID)"
