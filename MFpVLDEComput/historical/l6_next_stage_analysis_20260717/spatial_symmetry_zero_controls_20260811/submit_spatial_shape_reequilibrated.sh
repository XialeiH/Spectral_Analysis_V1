#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT=${1:?code root required}
SETUP_FILE=${2:?setup file required}
OUTPUT_ROOT=${3:?output root required}
mkdir -p "$OUTPUT_ROOT" "$OUTPUT_ROOT/logs"

ARRAY_JOB=$(sbatch --parsable \
  --account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --array=1-2%2 --cpus-per-task=8 --mem=48G --time=06:00:00 \
  --job-name=shape_reeq \
  --output="$OUTPUT_ROOT/logs/task_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/task_%A_%a.err" \
  --export=ALL,SPATIAL_CONTROL_SETUP="$SETUP_FILE",SPATIAL_CONTROL_OUTPUT="$OUTPUT_ROOT" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath(genpath('$CODE_ROOT')); run_spatial_shape_reequilibrated_task;\"")

PLOT_JOB=$(sbatch --parsable \
  --account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --dependency="afterok:$ARRAY_JOB" --cpus-per-task=2 --mem=12G --time=00:30:00 \
  --job-name=shape_reeq_plot \
  --output="$OUTPUT_ROOT/logs/plot_%j.out" \
  --error="$OUTPUT_ROOT/logs/plot_%j.err" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath(genpath('$CODE_ROOT')); plot_spatial_shape_reequilibrated_maps('$SETUP_FILE','$OUTPUT_ROOT');\"")

printf 'array_job\t%s\nplot_job\t%s\n' "$ARRAY_JOB" "$PLOT_JOB" | \
  tee "$OUTPUT_ROOT/submission_manifest.tsv"
