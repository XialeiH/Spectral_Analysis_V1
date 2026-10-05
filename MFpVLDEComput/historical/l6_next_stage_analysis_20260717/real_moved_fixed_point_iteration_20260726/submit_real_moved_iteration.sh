#!/bin/bash
set -euo pipefail

run_root=${1:?Usage: submit_real_moved_iteration.sh RUN_ROOT}
code_root="${run_root}/code"
setup_file="${run_root}/global_bifurcation_setup.mat"
output_root="${run_root}/results"
main_root=${REAL_MOVED_MAIN_ROOT:-/scratch/xh2906/NYU-Vision-2Drive-main}
h96_runtime=${REAL_MOVED_H96_RUNTIME:-/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle}

test -f "${setup_file}"
test -f "${code_root}/run_real_moved_iteration_task.m"
test -f "${code_root}/real_tuning_fixed_point.m"
test -f "${code_root}/real_tuning_true_jacobian.m"
test -f "${code_root}/l6ns_phi.m"
test -d "${main_root}"
test -d "${h96_runtime}"
mkdir -p "${output_root}/logs"

common_export="ALL,REAL_MOVED_CODE_ROOT=${code_root},REAL_MOVED_SETUP_FILE=${setup_file},REAL_MOVED_OUTPUT_ROOT=${output_root},REAL_MOVED_MAIN_ROOT=${main_root},REAL_MOVED_H96_RUNTIME=${h96_runtime}"
array_job=$(sbatch --parsable --array=1-10%10 \
    --output="${output_root}/logs/task_%A_%a.out" \
    --error="${output_root}/logs/task_%A_%a.err" \
    --export="${common_export}" \
    "${code_root}/run_real_moved_iteration_task.sbatch")
aggregate_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
    --output="${output_root}/logs/aggregate_%j.out" \
    --error="${output_root}/logs/aggregate_%j.err" \
    --export="${common_export}" \
    "${code_root}/aggregate_real_moved_iteration_results.sbatch")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\toutput_root\n' > "${manifest}"
printf 'iteration_array\t%s\t%s\n' "${array_job}" "${output_root}" >> "${manifest}"
printf 'aggregate\t%s\t%s\n' "${aggregate_job}" "${output_root}" >> "${manifest}"
printf 'array_job=%s\naggregate_job=%s\nmanifest=%s\n' \
    "${array_job}" "${aggregate_job}" "${manifest}"
