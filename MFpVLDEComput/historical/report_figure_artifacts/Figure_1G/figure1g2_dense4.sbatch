#!/bin/bash
#SBATCH --job-name=f1g2_d4
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-100%6
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
set -euo pipefail
module purge; module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g_dense_cold_equilibrium('$SOURCE_ROOT','$OUTPUT_ROOT',4,$SLURM_ARRAY_TASK_ID)"
