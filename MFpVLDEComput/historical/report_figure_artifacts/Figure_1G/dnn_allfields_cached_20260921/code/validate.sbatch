#!/bin/bash
#SBATCH --job-name=dnn_fields_pairs
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-220%16
set -euo pipefail
module purge
module load matlab/2025b
INDEX=$((SLURM_ARRAY_TASK_ID-1))
FIELDS=(6 8 10 20 30 40)
if (( INDEX < 100 )); then
    FIELD=4; REPEAT=$((INDEX+1))
else
    FIELD=${FIELDS[(INDEX-100)/20]}; REPEAT=$(((INDEX-100)%20+1))
fi
if (( REPEAT % 2 )); then
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','paged_combined','validation');"
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','reference','control');"
else
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','reference','control');"
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','paged_combined','validation');"
fi
