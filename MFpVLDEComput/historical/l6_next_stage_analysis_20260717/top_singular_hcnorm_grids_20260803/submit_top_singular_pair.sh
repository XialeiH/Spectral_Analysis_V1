#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_top_singular_pair.sh RUN_ROOT PAIR_TYPE}
pair_type=${2:?Usage: submit_top_singular_pair.sh RUN_ROOT PAIR_TYPE}
case "$pair_type" in
  l6_inhibition|l6_excitation|excitation_inhibition) ;;
  *) echo "Unsupported pair type: $pair_type" >&2; exit 1 ;;
esac

code_root="$run_root/code"
output_root="$run_root/$pair_type"
mkdir -p "$output_root/points" "$output_root/logs"
"$code_root/submit_top_singular_chunk.sh" "$run_root" "$pair_type" 1
