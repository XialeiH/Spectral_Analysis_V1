#!/bin/bash
set -euo pipefail
export OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1b_dnn_allfields_cached_20260921
export SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_field_scaling_20260825_001500
test ! -e "$OUTPUT_ROOT/submission.txt"
mkdir -p "$OUTPUT_ROOT"/{gate,validation,control,profile,logs}
GATE=$(sbatch --parsable --output="$OUTPUT_ROOT/logs/gate_%A_%a.log" "$OUTPUT_ROOT/code/gate.sbatch")
printf 'gate=%s\n' "$GATE" >> "$OUTPUT_ROOT/submission.txt"
VALIDATE=$(sbatch --parsable --dependency="afterok:$GATE" --output="$OUTPUT_ROOT/logs/validate_%A_%a.log" "$OUTPUT_ROOT/code/validate.sbatch")
printf 'validate=%s\n' "$VALIDATE" >> "$OUTPUT_ROOT/submission.txt"
cat "$OUTPUT_ROOT/submission.txt"
