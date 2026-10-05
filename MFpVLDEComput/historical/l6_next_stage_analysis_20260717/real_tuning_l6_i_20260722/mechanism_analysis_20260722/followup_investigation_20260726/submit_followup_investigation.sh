#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_followup_investigation.sh RUN_ROOT}
code_root="${run_root}/code"
output_root="${run_root}/results"
setup_file="${run_root}/global_bifurcation_setup.mat"
continuation_file="${run_root}/branch_continuation.mat"
main_root=${FOLLOWUP_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${FOLLOWUP_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}

test -f "${setup_file}"
test -f "${continuation_file}"
test -f "${code_root}/run_followup_branch_task.m"
test -d "${main_root}"
test -d "${h96_runtime}"
mkdir -p "${output_root}/logs"

common="ALL,FOLLOWUP_CODE_ROOT=${code_root},FOLLOWUP_SETUP_FILE=${setup_file},FOLLOWUP_CONTINUATION_FILE=${continuation_file},FOLLOWUP_OUTPUT_ROOT=${output_root},FOLLOWUP_MAIN_ROOT=${main_root},FOLLOWUP_H96_RUNTIME=${h96_runtime}"
branch=$(sbatch --parsable --array=1-28%8 --export="${common},FOLLOWUP_ENTRY=run_followup_branch_task" --output="${output_root}/logs/branch_%A_%a.out" --error="${output_root}/logs/branch_%A_%a.err" "${code_root}/followup_task.sbatch")
sensitivity=$(sbatch --parsable --array=1-9%8 --export="${common},FOLLOWUP_ENTRY=run_followup_sensitivity_task" --output="${output_root}/logs/sensitivity_%A_%a.out" --error="${output_root}/logs/sensitivity_%A_%a.err" "${code_root}/followup_task.sbatch")
interpolation=$(sbatch --parsable --array=1-3%3 --export="${common},FOLLOWUP_ENTRY=run_followup_interpolation_task" --output="${output_root}/logs/interpolation_%A_%a.out" --error="${output_root}/logs/interpolation_%A_%a.err" "${code_root}/followup_task.sbatch")
threshold=$(sbatch --parsable --array=1-12%8 --export="${common},FOLLOWUP_ENTRY=run_followup_threshold_task" --output="${output_root}/logs/threshold_%A_%a.out" --error="${output_root}/logs/threshold_%A_%a.err" "${code_root}/followup_task.sbatch")
ode=$(sbatch --parsable --export="${common},FOLLOWUP_ENTRY=run_followup_ode_hidden_mode" --output="${output_root}/logs/ode_%j.out" --error="${output_root}/logs/ode_%j.err" "${code_root}/followup_task.sbatch")
dependency="afterok:${branch}:${sensitivity}:${interpolation}:${threshold}:${ode}"
aggregate=$(sbatch --parsable --dependency="${dependency}" --export="${common},FOLLOWUP_ENTRY=aggregate_followup_results" --output="${output_root}/logs/aggregate_%j.out" --error="${output_root}/logs/aggregate_%j.err" "${code_root}/followup_task.sbatch")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\toutput_root\n' > "${manifest}"
for item in "branch:${branch}" "sensitivity:${sensitivity}" "interpolation:${interpolation}" "threshold:${threshold}" "ode:${ode}" "aggregate:${aggregate}"; do
    printf '%s\t%s\t%s\n' "${item%%:*}" "${item##*:}" "${output_root}" >> "${manifest}"
done
cat "${manifest}"
