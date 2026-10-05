#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_pair_grid_chunk.sh RUN_ROOT PAIR_TYPE SETUP_FILE START_TASK}
pair_type=${2:?Usage: submit_pair_grid_chunk.sh RUN_ROOT PAIR_TYPE SETUP_FILE START_TASK}
setup_file=${3:?Usage: submit_pair_grid_chunk.sh RUN_ROOT PAIR_TYPE SETUP_FILE START_TASK}
start_task=${4:?Usage: submit_pair_grid_chunk.sh RUN_ROOT PAIR_TYPE SETUP_FILE START_TASK}
total_count=6561
chunk_size=2800
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
main_root=${PAIR_GRID_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${PAIR_GRID_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}
common="ALL,PAIR_GRID_CODE_ROOT=$code_root,PAIR_GRID_SETUP_FILE=$setup_file,PAIR_GRID_OUTPUT_ROOT=$output_root,PAIR_GRID_MAIN_ROOT=$main_root,PAIR_GRID_H96_RUNTIME=$h96_runtime,PAIR_GRID_PAIR=$pair_type"

array_job=$(sbatch --parsable --array="1-${count}%24" \
  --partition=cs,cpu_short \
  --output=/dev/null --error=/dev/null \
  --export="$common,PAIR_GRID_TASK_OFFSET=$offset" \
  "$code_root/run_pair_grid_task.sbatch")

next_task=$((end_task+1))
if (( next_task <= total_count )); then
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/controller_%j.out" \
    --error="$output_root/logs/controller_%j.err" \
    --export="ALL,PAIR_GRID_CODE_ROOT=$code_root,PAIR_GRID_RUN_ROOT=$run_root,PAIR_GRID_PAIR=$pair_type,PAIR_GRID_SETUP_FILE=$setup_file,PAIR_GRID_NEXT_START=$next_task" \
    "$code_root/continue_pair_grid.sbatch")
  followup_stage="controller"
else
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/aggregate_%j.out" \
    --error="$output_root/logs/aggregate_%j.err" \
    --export="$common" "$code_root/aggregate_pair_grid.sbatch")
  followup_stage="aggregate"
fi

manifest="$run_root/submission_manifest.tsv"
{
  flock 9
  printf 'full_grid_array\t%d\t%d\t%s\t%s\t%s\n' \
    "$start_task" "$end_task" "$array_job" "$pair_type" "$output_root" >> "$manifest"
  printf '%s\t%d\t%d\t%s\t%s\t%s\n' \
    "$followup_stage" "$next_task" "$next_task" "$followup_job" "$pair_type" "$output_root" >> "$manifest"
} 9>"$run_root/submission_manifest.lock"
printf 'array_job=%s range=%d-%d %s_job=%s pair=%s\n' \
  "$array_job" "$start_task" "$end_task" "$followup_stage" "$followup_job" "$pair_type"
