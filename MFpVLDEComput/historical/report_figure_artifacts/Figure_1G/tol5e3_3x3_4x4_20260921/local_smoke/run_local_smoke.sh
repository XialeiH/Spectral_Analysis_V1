#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SOURCE_ROOT="$HERE/source"
OUTPUT_ROOT="$HERE/output"
# Preserve the benchmark implementation, using this Mac's 12 logical CPUs.
export SLURM_CPUS_PER_TASK=12
for REPEAT in 1 2 3; do
    for FIELD in 3 4; do
        if (( REPEAT % 2 == 1 )); then
            METHODS=("CG" "DNN surrogate")
        else
            METHODS=("DNN surrogate" "CG")
        fi
        for METHOD in "${METHODS[@]}"; do
            KEY="${METHOD// /_}"
            LOG="$OUTPUT_ROOT/logs/${KEY}_${FIELD}_${REPEAT}.log"
            printf 'Local smoke: field=%s repeat=%s method=%s\n' "$FIELD" "$REPEAT" "$METHOD"
            if ! /Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g2_tol5e3('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','$METHOD')" > "$LOG" 2>&1; then
                cat "$LOG"
                exit 1
            fi
            tail -n 3 "$LOG"
        done
    done
done

