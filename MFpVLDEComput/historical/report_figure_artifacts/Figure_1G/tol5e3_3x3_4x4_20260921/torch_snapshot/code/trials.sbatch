#!/bin/bash
#SBATCH --job-name=f1b_tol5e3
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-400%3
set -euo pipefail
module purge
module load matlab/2025b
INDEX=$((SLURM_ARRAY_TASK_ID-1))
FIELD=$((3+INDEX/200))
REPEAT=$((INDEX%100+1))
if (( INDEX%200 < 100 )); then
    METHOD='CG'
else
    METHOD='DNN surrogate'
fi
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g2_tol5e3('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','$METHOD')"

