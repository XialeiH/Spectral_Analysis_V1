#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_real_tuning.sh RUN_ROOT [smoke]}
mode=${2:-full}
code_root="$run_root/code"
output_root="$run_root/results"
setup_file="$run_root/global_bifurcation_setup.mat"
main_root=${REAL_TUNING_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${REAL_TUNING_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}

if [[ "$mode" == smoke ]]; then
  smoke=1
  count=18
else
  smoke=0
  count=82
fi
mkdir -p "$output_root/logs"
common="ALL,REAL_TUNING_CODE_ROOT=${code_root},REAL_TUNING_SETUP_FILE=${setup_file},REAL_TUNING_OUTPUT_ROOT=${output_root},REAL_TUNING_MAIN_ROOT=${main_root},REAL_TUNING_H96_RUNTIME=${h96_runtime},REAL_TUNING_SMOKE=${smoke}"
array_job=$(sbatch --parsable --array="1-${count}%48" \
  --output="$output_root/logs/task_%A_%a.out" \
  --error="$output_root/logs/task_%A_%a.err" \
  --export="$common" "$code_root/run_real_tuning_task.sbatch")
aggregate_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
  --output="$output_root/logs/aggregate_%j.out" \
  --error="$output_root/logs/aggregate_%j.err" \
  --export="$common" "$code_root/aggregate_real_tuning_results.sbatch")
manifest="$run_root/submission_manifest_${mode}.tsv"
printf 'stage\tjob_id\toutput_root\n' > "$manifest"
printf 'beta_array\t%s\t%s\naggregate\t%s\t%s\n' \
  "$array_job" "$output_root" "$aggregate_job" "$output_root" >> "$manifest"
printf 'array_job=%s\naggregate_job=%s\nmanifest=%s\n' \
  "$array_job" "$aggregate_job" "$manifest"
