#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_global_offset.sh RUN_ROOT [smoke]}
mode=${2:-full}
code_root="$run_root/code"
output_root="$run_root/results"
source_root=${GLOBAL_OFFSET_SOURCE_ROOT:-/scratch/xh2906/librarySCI_runs/l6_i_offset_full_20260722_192837}
fpp_root=${GLOBAL_OFFSET_FPP_ROOT:-/scratch/xh2906/librarySCI_runs/fpp_offset_comparison_full_correctJI_20260722_200304/results}
main_root=${GLOBAL_OFFSET_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${GLOBAL_OFFSET_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}
if [[ "$mode" == smoke ]]; then smoke=1; count=3; throttle=3; else smoke=0; count=126; throttle=40; fi
mkdir -p "$output_root/logs"
common="ALL,GLOBAL_OFFSET_CODE_ROOT=$code_root,GLOBAL_OFFSET_SETUP_FILE=$run_root/global_bifurcation_setup.mat,GLOBAL_OFFSET_OUTPUT_ROOT=$output_root,GLOBAL_OFFSET_OLD_SUMMARY=$source_root/results/offset_summary.tsv,GLOBAL_OFFSET_FPP_SUMMARY=$fpp_root/fpp_offset_summary.tsv,GLOBAL_OFFSET_MAIN_ROOT=$main_root,GLOBAL_OFFSET_H96_RUNTIME=$h96_runtime,GLOBAL_OFFSET_SMOKE=$smoke"
array_job=$(sbatch --parsable --array="1-${count}%${throttle}" \
  --output="$output_root/logs/task_%A_%a.out" \
  --error="$output_root/logs/task_%A_%a.err" \
  --export="$common" "$code_root/run_global_offset_task.sbatch")
aggregate_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
  --output="$output_root/logs/aggregate_%j.out" \
  --error="$output_root/logs/aggregate_%j.err" \
  --export="$common" "$code_root/aggregate_global_offset.sbatch")
manifest="$run_root/submission_manifest_${mode}.tsv"
printf 'stage\tjob_id\toutput_root\n' > "$manifest"
printf 'global_offset_array\t%s\t%s\naggregate\t%s\t%s\n' \
  "$array_job" "$output_root" "$aggregate_job" "$output_root" >> "$manifest"
printf 'array_job=%s\naggregate_job=%s\nmanifest=%s\n' \
  "$array_job" "$aggregate_job" "$manifest"
