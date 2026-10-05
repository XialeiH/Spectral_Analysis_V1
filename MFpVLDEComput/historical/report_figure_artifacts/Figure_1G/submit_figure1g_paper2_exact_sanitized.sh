#!/bin/bash
set -euo pipefail

BASE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829
EXACT_ROOT=$BASE_ROOT/snn_paper2_original_exact_20260901_160003

test -f "$EXACT_ROOT/code/figure1g_paper2_exact_plot.sbatch"
mkdir -p "$EXACT_ROOT/logs" "$EXACT_ROOT/results"

ARRAY_ID=$(sbatch --parsable \
    --export=ALL,BASE_ROOT="$BASE_ROOT",EXACT_ROOT="$EXACT_ROOT" \
    --output="$EXACT_ROOT/logs/sanitized_%A_%a.out" \
    --error="$EXACT_ROOT/logs/sanitized_%A_%a.err" \
    "$EXACT_ROOT/code/figure1g_paper2_exact_trials_sanitized.sbatch")
PLOT_ID=$(sbatch --parsable \
    --dependency=afterok:"$ARRAY_ID" \
    --export=ALL,BASE_ROOT="$BASE_ROOT",EXACT_ROOT="$EXACT_ROOT" \
    --output="$EXACT_ROOT/logs/sanitized_plot_%j.out" \
    --error="$EXACT_ROOT/logs/sanitized_plot_%j.err" \
    "$EXACT_ROOT/code/figure1g_paper2_exact_plot.sbatch")

printf 'sanitized_retry_array\t%s\t%s\n' "$ARRAY_ID" "$(date -Iseconds)" >> "$EXACT_ROOT/submission_manifest.tsv"
printf 'sanitized_retry_plot\t%s\t%s\n' "$PLOT_ID" "$(date -Iseconds)" >> "$EXACT_ROOT/submission_manifest.tsv"
printf 'ARRAY_ID=%s\nPLOT_ID=%s\n' "$ARRAY_ID" "$PLOT_ID"
squeue -j "$ARRAY_ID,$PLOT_ID" -o '%.18i %.9P %.28j %.2t %.10M %.6D %R'
