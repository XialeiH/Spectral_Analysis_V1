#!/bin/bash
#SBATCH --job-name=f1g2_method
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
module purge; module load matlab/2025b
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
FIELD=$(printf '%s\n' "$LINE" | cut -f1)
REPEAT=$(printf '%s\n' "$LINE" | cut -f2)
METHOD=$(printf '%s\n' "$LINE" | cut -f3)
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g2_method('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','$METHOD')"
