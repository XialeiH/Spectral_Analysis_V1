#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_l6_ei_grid.sh RUN_ROOT}
mkdir -p "$run_root/results/points" "$run_root/results/logs"
manifest="$run_root/submission_manifest.tsv"
printf 'stage\tstart_task\tend_task\tjob_id\toutput_root\n' > "$manifest"
"$run_root/code/submit_l6_ei_grid_chunk.sh" "$run_root" 1
printf 'manifest=%s\n' "$manifest"
