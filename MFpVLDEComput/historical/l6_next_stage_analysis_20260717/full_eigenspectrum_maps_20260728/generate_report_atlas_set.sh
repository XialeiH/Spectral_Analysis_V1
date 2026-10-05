#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 SOURCE_ROOT OUTPUT_ROOT [START_POPULATION:JACOBIAN]" >&2
    exit 2
fi

source_root="$1"
output_root="$2"
start_key="${3:-}"
matlab_bin="/Applications/MATLAB_R2025b.app/bin/matlab"
started=0
if [[ -z "$start_key" ]]; then
    started=1
fi

for jacobian in J_full J_6 J_I; do
    for population in E S C I; do
        current_key="${population}:${jacobian}"
        if (( started == 0 )); then
            if [[ "$current_key" != "$start_key" ]]; then
                continue
            fi
            started=1
        fi

        set +e
        "$matlab_bin" -batch \
            "replot_report_full_eigenspectrum_atlases('$source_root','$output_root','$population','$jacobian');"
        matlab_status=$?
        set -e

        pdf_file="$output_root/${population}_population_${jacobian}_16condition_full_eigenspectrum_maps.pdf"
        pdf_size=0
        if [[ -f "$pdf_file" ]]; then
            pdf_size=$(stat -f '%z' "$pdf_file")
        fi
        if (( pdf_size < 100000 )); then
            echo "Missing or truncated output after $current_key" >&2
            exit 1
        fi
        if (( matlab_status != 0 )); then
            echo "MATLAB exited after saving validated output for $current_key; continuing."
        fi
    done
done
