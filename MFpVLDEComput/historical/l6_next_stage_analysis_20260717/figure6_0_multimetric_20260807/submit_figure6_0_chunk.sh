#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_figure6_0_chunk.sh RUN_ROOT START_TASK}
start_task=${2:?Usage: submit_figure6_0_chunk.sh RUN_ROOT START_TASK}
total_count=10404
chunk_size=2600
if (( start_task < 1 || start_task > total_count )); then
  echo "Invalid start task: $start_task" >&2
  exit 1
fi
remaining=$((total_count-start_task+1))
if (( remaining < chunk_size )); then count=$remaining; else count=$chunk_size; fi
end_task=$((start_task+count-1))
offset=$((start_task-1))

code_root="$run_root/code"
output_root="$run_root/results"
main_root=/scratch/xh2906/NYU-Vision-2Drive-main
h96_runtime=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle_extended_corrected_h96_ikp7_l6p5_matchedlib_20260702_010052
common="ALL,FIG60_RUN_ROOT=$run_root,FIG60_CODE_ROOT=$code_root,FIG60_SETUP_ROOT=$run_root/setup,FIG60_OUTPUT_ROOT=$output_root,FIG60_FIXED_TOLERANCE=1e-4,FIG60_MAIN_ROOT=$main_root,FIG60_H96_RUNTIME=$h96_runtime"

array_job=$(sbatch --parsable --array="1-${count}%48" \
  --partition=cs,cpu_short --output=/dev/null --error=/dev/null \
  --export="$common,FIG60_TASK_OFFSET=$offset" \
  "$code_root/run_figure6_0_grid_task.sbatch")
next_task=$((end_task+1))
if (( next_task <= total_count )); then
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/controller_%j.out" \
    --error="$output_root/logs/controller_%j.err" \
    --export="ALL,FIG60_RUN_ROOT=$run_root,FIG60_CODE_ROOT=$code_root,FIG60_NEXT_START=$next_task" \
    "$code_root/continue_figure6_0_chunk.sbatch")
  followup_stage=controller
else
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/finalize_%j.out" \
    --error="$output_root/logs/finalize_%j.err" \
    --export="$common" "$code_root/finalize_figure6_0.sbatch")
  followup_stage=finalize
fi
{
  flock 9
  printf 'array\t%d\t%d\t%s\n' "$start_task" "$end_task" "$array_job" >&9
  printf '%s\t%d\t%s\n' "$followup_stage" "$next_task" "$followup_job" >&9
} 9>>"$run_root/submission_manifest.tsv"
printf 'array_job=%s range=%d-%d %s_job=%s\n' \
  "$array_job" "$start_task" "$end_task" "$followup_stage" "$followup_job"
