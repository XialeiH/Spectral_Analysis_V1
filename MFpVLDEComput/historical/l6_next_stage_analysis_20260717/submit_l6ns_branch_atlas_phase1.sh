#!/usr/bin/env bash
set -euo pipefail

root=${1:?Usage: submit_l6ns_branch_atlas_phase1.sh RUN_ROOT}
catalog="$root/branch_atlas_20260803/atlas_seed_catalog.tsv"
mkdir -p "$root/logs" "$root/branch_atlas_20260803/continuation_tasks"
seed_count=$(awk -F '\t' 'NR>1 && ($5=="census_root" || $5=="branch_intersection" || $5=="nearsilent_family" || $5=="nearsilent_fold_seed") {n++} END {print n+0}' "$catalog")
task_count=$((2*seed_count))
if (( task_count == 0 )); then
  echo "No eligible seeds in $catalog" >&2
  exit 1
fi
export L6NS_ATLAS_ROOT="$root"
job_id=$(sbatch --parsable --array="1-${task_count}%16" --export=ALL,L6NS_ATLAS_ROOT="$root" "$root/run_l6ns_branch_atlas_task.sbatch")
aggregate_job_id=$(sbatch --parsable --dependency="afterok:${job_id}" --export=ALL,L6NS_ATLAS_ROOT="$root" "$root/aggregate_l6ns_branch_atlas_phase1.sbatch")
manifest="$root/branch_atlas_20260803/phase1_submission_manifest.tsv"
printf 'job_id\taggregate_job_id\tseed_count\ttask_count\trun_root\n%s\t%s\t%d\t%d\t%s\n' "$job_id" "$aggregate_job_id" "$seed_count" "$task_count" "$root" > "$manifest"
printf 'JOB_ID=%s\nAGGREGATE_JOB_ID=%s\nSEEDS=%d\nTASKS=%d\nMANIFEST=%s\n' "$job_id" "$aggregate_job_id" "$seed_count" "$task_count" "$manifest"
