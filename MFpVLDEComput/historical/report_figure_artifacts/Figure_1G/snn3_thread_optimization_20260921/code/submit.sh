#!/bin/bash
set -euo pipefail
export OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1b_snn3_thread_optimization_20260921
test ! -e "$OUTPUT_ROOT/submission.txt"
mkdir -p "$OUTPUT_ROOT/results" "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/trials"
echo "b3958be55c41dba44f674de73e3dfb8bf37e2b552c7f4b145b824872197f4c6b  $OUTPUT_ROOT/code/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
bash -n "$OUTPUT_ROOT/code/study.sbatch"
probe=$(sbatch --parsable --array=1-15%3 --output="$OUTPUT_ROOT/logs/probe_%A_%a.out" --error="$OUTPUT_ROOT/logs/probe_%A_%a.err" "$OUTPUT_ROOT/code/study.sbatch" probe)
printf 'probe=%s\n' "$probe" >> "$OUTPUT_ROOT/submission.txt"
selection=$(sbatch --parsable --dependency="afterok:$probe" --output="$OUTPUT_ROOT/logs/select_%j.out" --error="$OUTPUT_ROOT/logs/select_%j.err" "$OUTPUT_ROOT/code/study.sbatch" select)
printf 'selection=%s\n' "$selection" >> "$OUTPUT_ROOT/submission.txt"
validation=$(sbatch --parsable --array=1-5%3 --dependency="afterok:$selection" --output="$OUTPUT_ROOT/logs/validation_%A_%a.out" --error="$OUTPUT_ROOT/logs/validation_%A_%a.err" "$OUTPUT_ROOT/code/study.sbatch" validate)
printf 'validation=%s\n' "$validation" >> "$OUTPUT_ROOT/submission.txt"
cat "$OUTPUT_ROOT/submission.txt"
