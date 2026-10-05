#!/bin/bash
set -euo pipefail
export OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1b_dnn3_impl_optimization_20260921
export SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_field_scaling_20260825_001500
mkdir -p "$OUTPUT_ROOT"/{probe,profile,validation,control,logs}
if [[ -s "$OUTPUT_ROOT/submission.txt" ]]; then
    printf 'Already submitted; not submitting duplicates.\n'
    cat "$OUTPUT_ROOT/submission.txt"
    exit 0
fi
PROBE=$(sbatch --parsable --output="$OUTPUT_ROOT/logs/probe_%A_%a.log" "$OUTPUT_ROOT/code/probe.sbatch")
printf 'probe=%s\n' "$PROBE" >> "$OUTPUT_ROOT/submission.txt"
SELECT=$(sbatch --parsable --dependency="afterok:$PROBE" --output="$OUTPUT_ROOT/logs/select_%j.log" "$OUTPUT_ROOT/code/select.sbatch")
printf 'select=%s\n' "$SELECT" >> "$OUTPUT_ROOT/submission.txt"
VALIDATE=$(sbatch --parsable --dependency="afterok:$SELECT" --output="$OUTPUT_ROOT/logs/validate_%A_%a.log" "$OUTPUT_ROOT/code/validate.sbatch")
printf 'validate=%s\n' "$VALIDATE" >> "$OUTPUT_ROOT/submission.txt"
cat "$OUTPUT_ROOT/submission.txt"
