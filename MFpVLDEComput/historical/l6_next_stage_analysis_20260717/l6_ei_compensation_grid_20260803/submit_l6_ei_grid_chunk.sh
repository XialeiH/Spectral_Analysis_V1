#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_l6_ei_grid_chunk.sh RUN_ROOT START_TASK}
start_task=${2:?Usage: submit_l6_ei_grid_chunk.sh RUN_ROOT START_TASK}
if [[ -f "$run_root/grid_config.sh" ]]; then
  source "$run_root/grid_config.sh"
fi
: "${L6EI_GRID_STEP:=0.005}"
: "${L6EI_BETA6_MAXIMUM:=0.35}"
: "${L6EI_BETAEI_MAXIMUM:=0.2}"
: "${L6EI_TOTAL_GRID_COUNT:=2911}"
: "${L6EI_CHUNK_SIZE:=2600}"
: "${L6EI_TASK_LIST_FILE:=$run_root/missing_task_ids.tsv}"
if [[ -f "$L6EI_TASK_LIST_FILE" ]]; then
  total_count=$(wc -l < "$L6EI_TASK_LIST_FILE")
else
  total_count=$L6EI_TOTAL_GRID_COUNT
  L6EI_TASK_LIST_FILE=
fi
chunk_size=$L6EI_CHUNK_SIZE
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
h96_runtime=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle
common="ALL,L6EI_CODE_ROOT=$code_root,L6EI_SETUP_FILE=$run_root/global_bifurcation_setup.mat,L6EI_OUTPUT_ROOT=$output_root,L6EI_MAIN_ROOT=$main_root,L6EI_H96_RUNTIME=$h96_runtime,L6EI_TASK_LIST_FILE=$L6EI_TASK_LIST_FILE,L6EI_GRID_STEP=$L6EI_GRID_STEP,L6EI_BETA6_MAXIMUM=$L6EI_BETA6_MAXIMUM,L6EI_BETAEI_MAXIMUM=$L6EI_BETAEI_MAXIMUM,L6EI_TOTAL_GRID_COUNT=$L6EI_TOTAL_GRID_COUNT"

array_job=$(sbatch --parsable --partition=cs,cpu_short --array="1-${count}%48" \
  --output=/dev/null --error=/dev/null \
  --export="$common,L6EI_TASK_OFFSET=$offset" \
  "$code_root/run_l6_ei_grid_task.sbatch")

next_task=$((end_task+1))
if (( next_task <= total_count )); then
  followup_job=$(sbatch --parsable --partition=cs,cpu_short \
    --dependency="afterok:${array_job}" \
    --output="$output_root/logs/controller_%j.out" \
    --error="$output_root/logs/controller_%j.err" \
    --export="ALL,L6EI_CODE_ROOT=$code_root,L6EI_RUN_ROOT=$run_root,L6EI_NEXT_START=$next_task" \
    "$code_root/continue_l6_ei_grid.sbatch")
  followup_stage="controller"
else
  followup_job=$(sbatch --parsable --partition=cs,cpu_short \
    --dependency="afterok:${array_job}" \
    --output="$output_root/logs/aggregate_%j.out" \
    --error="$output_root/logs/aggregate_%j.err" \
    --export="$common" "$code_root/aggregate_l6_ei_grid.sbatch")
  followup_stage="aggregate"
fi

manifest="$run_root/submission_manifest.tsv"
{
  flock 9
  printf 'grid_array\t%d\t%d\t%s\t%s\n' \
    "$start_task" "$end_task" "$array_job" "$output_root" >> "$manifest"
  printf '%s\t%d\t%d\t%s\t%s\n' \
    "$followup_stage" "$next_task" "$next_task" "$followup_job" "$output_root" >> "$manifest"
} 9>"$run_root/submission_manifest.lock"
printf 'array_job=%s range=%d-%d %s_job=%s\n' \
  "$array_job" "$start_task" "$end_task" "$followup_stage" "$followup_job"
