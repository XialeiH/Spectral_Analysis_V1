#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_pair_grid.sh RUN_ROOT PAIR_TYPE SETUP_FILE}
pair_type=${2:?Usage: submit_pair_grid.sh RUN_ROOT PAIR_TYPE SETUP_FILE}
setup_file=${3:?Usage: submit_pair_grid.sh RUN_ROOT PAIR_TYPE SETUP_FILE}
if [[ "$pair_type" != l6_excitation && "$pair_type" != excitation_inhibition ]]; then
  echo "Unsupported pair type: $pair_type" >&2
  exit 1
fi
mkdir -p "$run_root/results/points" "$run_root/results/logs"
manifest="$run_root/submission_manifest.tsv"
printf 'stage\tstart_task\tend_task\tjob_id\tpair_type\toutput_root\n' > "$manifest"
"$run_root/code/submit_pair_grid_chunk.sh" "$run_root" "$pair_type" "$setup_file" 1
printf 'manifest=%s\n' "$manifest"
