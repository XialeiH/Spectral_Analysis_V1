#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_top_singular_chunk.sh RUN_ROOT PAIR_TYPE START_TASK}
pair_type=${2:?Usage: submit_top_singular_chunk.sh RUN_ROOT PAIR_TYPE START_TASK}
start_task=${3:?Usage: submit_top_singular_chunk.sh RUN_ROOT PAIR_TYPE START_TASK}
total_count=2601
chunk_size=400
if (( start_task < 1 || start_task > total_count )); then
  echo "Invalid start task: $start_task" >&2
  exit 1
fi
remaining=$((total_count-start_task+1))
if (( remaining < chunk_size )); then count=$remaining; else count=$chunk_size; fi
end_task=$((start_task+count-1))
offset=$((start_task-1))

code_root="$run_root/code"
output_root="$run_root/$pair_type"
setup_file=/scratch/xh2906/librarySCI_runs/l6_i_canonicalHC_full_grid_20260728_001/global_bifurcation_setup.mat
baseline_mode_file="$run_root/baseline/baseline_top_singular_mode.mat"
main_root=/scratch/xh2906/NYU-Vision-2Drive-main
h96_runtime=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle
common="ALL,TSHC_RUN_ROOT=$run_root,TSHC_CODE_ROOT=$code_root,TSHC_SETUP_FILE=$setup_file,TSHC_BASELINE_MODE_FILE=$baseline_mode_file,TSHC_OUTPUT_ROOT=$output_root,TSHC_PAIR=$pair_type,TSHC_FIXED_TOLERANCE=1e-4,TSHC_MAIN_ROOT=$main_root,TSHC_H96_RUNTIME=$h96_runtime"

array_job=$(sbatch --parsable --array="1-${count}%24" \
  --partition=cs,cpu_short --output=/dev/null --error=/dev/null \
  --export="$common,TSHC_TASK_OFFSET=$offset" \
  "$code_root/run_top_singular_grid_task.sbatch")
next_task=$((end_task+1))
if (( next_task <= total_count )); then
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/controller_%j.out" \
    --error="$output_root/logs/controller_%j.err" \
    --export="ALL,TSHC_RUN_ROOT=$run_root,TSHC_CODE_ROOT=$code_root,TSHC_PAIR=$pair_type,TSHC_NEXT_START=$next_task" \
    "$code_root/continue_top_singular_chunk.sbatch")
  followup_stage=controller
else
  followup_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="$output_root/logs/aggregate_%j.out" \
    --error="$output_root/logs/aggregate_%j.err" \
    --export="$common" "$code_root/aggregate_top_singular_pair.sbatch")
  followup_stage=aggregate
fi
{
  flock 9
  printf 'array\t%s\t%d\t%d\t%s\n' \
    "$pair_type" "$start_task" "$end_task" "$array_job" >&9
  printf '%s\t%s\t%d\t%d\t%s\n' \
    "$followup_stage" "$pair_type" "$next_task" "$next_task" "$followup_job" >&9
} 9>>"$run_root/submission_manifest.tsv"
printf 'pair=%s array_job=%s range=%d-%d %s_job=%s\n' \
  "$pair_type" "$array_job" "$start_task" "$end_task" \
  "$followup_stage" "$followup_job"
