#!/bin/bash
#SBATCH --job-name=f1b_tol5e3_large
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-210%16
set -euo pipefail
module purge
module load matlab/2025b
INDEX=$((SLURM_ARRAY_TASK_ID-1))
FIELDS=(6 8 10 20 30 40)
if (( INDEX < 80 )); then
    FIELD=${FIELDS[INDEX/20]}
    REPEAT=$((INDEX%20+1))
    METHOD='CG'
elif (( INDEX < 90 )); then
    FIELD=${FIELDS[4+(INDEX-80)/5]}
    REPEAT=$(((INDEX-80)%5+1))
    METHOD='CG'
else
    INDEX=$((INDEX-90))
    FIELD=${FIELDS[INDEX/20]}
    REPEAT=$((INDEX%20+1))
    METHOD='DNN surrogate'
fi
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g2_tol5e3('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','$METHOD')"
