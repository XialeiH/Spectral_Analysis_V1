#!/bin/bash
#SBATCH --job-name=f1g2_scale
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-30%6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
module purge; module load matlab/2025b
FIELDS=(6 8 10 20 30 40)
FIELD=${FIELDS[$(( (SLURM_ARRAY_TASK_ID-1) % 6 ))]}
REPEAT=$(( (SLURM_ARRAY_TASK_ID-1) / 6 + 1 ))
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g_dense_cold_equilibrium('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT')"
