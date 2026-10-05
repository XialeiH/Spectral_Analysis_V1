#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_fpp_offset_comparison.sh RUN_ROOT [smoke]}
mode=${2:-full}
code_root="$run_root/code"
output_root="$run_root/results"
setup_file="$run_root/global_bifurcation_setup.mat"
theory_file="$run_root/opposite_effect_summary.tsv"
real_root=${FPP_OFFSET_REAL_ROOT:-/scratch/xh2906/librarySCI_runs/l6_i_offset_full_20260722_192837/results}
main_root=${FPP_OFFSET_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${FPP_OFFSET_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}
if [[ "$mode" == smoke ]]; then smoke=1; count=3; else smoke=0; count=21; fi
mkdir -p "$output_root/logs"
common="ALL,FPP_OFFSET_CODE_ROOT=$code_root,FPP_OFFSET_SETUP_FILE=$setup_file,FPP_OFFSET_THEORY_FILE=$theory_file,FPP_OFFSET_OUTPUT_ROOT=$output_root,FPP_OFFSET_REAL_SUMMARY=$real_root/offset_summary.tsv,FPP_OFFSET_REAL_REGRESSION=$real_root/offset_regression.tsv,FPP_OFFSET_MAIN_ROOT=$main_root,FPP_OFFSET_H96_RUNTIME=$h96_runtime,FPP_OFFSET_SMOKE=$smoke"
array_job=$(sbatch --parsable --array="1-${count}%21" \
  --output="$output_root/logs/task_%A_%a.out" \
  --error="$output_root/logs/task_%A_%a.err" \
  --export="$common" "$code_root/run_fpp_offset_task.sbatch")
aggregate_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
  --output="$output_root/logs/aggregate_%j.out" \
  --error="$output_root/logs/aggregate_%j.err" \
  --export="$common" "$code_root/aggregate_fpp_offset_comparison.sbatch")
manifest="$run_root/submission_manifest_${mode}.tsv"
printf 'stage\tjob_id\toutput_root\n' > "$manifest"
printf 'fpp_offset_array\t%s\t%s\naggregate\t%s\t%s\n' \
  "$array_job" "$output_root" "$aggregate_job" "$output_root" >> "$manifest"
printf 'array_job=%s\naggregate_job=%s\nmanifest=%s\n' \
  "$array_job" "$aggregate_job" "$manifest"
