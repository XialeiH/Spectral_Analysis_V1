#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT=${1:?code root required}
SETUP_FILE=${2:?setup file required}
OUTPUT_ROOT=${3:?output root required}
mkdir -p "$OUTPUT_ROOT" "$OUTPUT_ROOT/logs"

ARRAY_JOB=$(sbatch --parsable \
  --account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --array=1-17%8 --cpus-per-task=8 --mem=48G --time=04:00:00 \
  --job-name=fine_alpha_overlap \
  --output="$OUTPUT_ROOT/logs/task_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/task_%A_%a.err" \
  --export=ALL,SPATIAL_CONTROL_SETUP="$SETUP_FILE",SPATIAL_CONTROL_OUTPUT="$OUTPUT_ROOT" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath(genpath('$CODE_ROOT')); run_fine_alpha_overlap_task;\"")

AGG_JOB=$(sbatch --parsable \
  --account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --dependency="afterok:$ARRAY_JOB" --cpus-per-task=4 --mem=24G --time=00:30:00 \
  --job-name=fine_alpha_agg \
  --output="$OUTPUT_ROOT/logs/aggregate_%j.out" \
  --error="$OUTPUT_ROOT/logs/aggregate_%j.err" \
  --export=ALL,SPATIAL_CONTROL_OUTPUT="$OUTPUT_ROOT" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath(genpath('$CODE_ROOT')); aggregate_fine_alpha_overlap;\"")

printf 'array_job\t%s\naggregate_job\t%s\n' "$ARRAY_JOB" "$AGG_JOB" | tee "$OUTPUT_ROOT/submission_manifest.tsv"
