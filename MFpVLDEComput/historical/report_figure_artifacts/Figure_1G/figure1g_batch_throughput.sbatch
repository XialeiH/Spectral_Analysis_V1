#!/bin/bash
#SBATCH --job-name=f1g_batch
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-12%6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/figure1g_evaluator_suite_20260829/logs/batch_%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/figure1g_evaluator_suite_20260829/logs/batch_%A_%a.err
set -euo pipefail
SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_field_scaling_20260825_001500
OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_evaluator_suite_20260829
FIELDS=(4 10 20 40); FIELD=${FIELDS[$(( (SLURM_ARRAY_TASK_ID-1) % 4 ))]}; REPEAT=$(( (SLURM_ARRAY_TASK_ID-1) / 4 + 1 ))
mkdir -p "$OUTPUT_ROOT/batch_throughput" "$OUTPUT_ROOT/logs"
module purge; module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g_batch_throughput('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT')"
