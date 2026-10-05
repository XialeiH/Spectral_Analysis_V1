#!/bin/bash
#SBATCH --job-name=dnn_fields_gate
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-14%7
set -euo pipefail
module purge
module load matlab/2025b
FIELDS=(4 6 8 10 20 30 40)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
FIELD=${FIELDS[INDEX/2]}
REPEAT=$((1+3*(INDEX%2)))
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','paged_combined','gate');"
