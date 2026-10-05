#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 RUN_ROOT EXISTING_DATA_ROOT LOCAL_FIGURE_DESTINATION [PROFILE]" >&2
    exit 2
fi

run_root=$1
existing_data_root=$2
local_figure_destination=$3
profile=${4:-initial}
case "${profile}" in
    initial) task_count=11 ;;
    remaining) task_count=8 ;;
    *) echo "Unsupported profile: ${profile}" >&2; exit 2 ;;
esac
code_root="${run_root}/code"
data_root="${run_root}/data"
figure_root="${run_root}/figures"
log_root="${run_root}/logs"
mkdir -p "${data_root}" "${figure_root}" "${log_root}"

export CLUSTER_CODE_ROOT="${code_root}"
export CLUSTER_EXISTING_DATA_ROOT="${existing_data_root}"
export CLUSTER_OUTPUT_ROOT="${data_root}"
export CLUSTER_PROFILE="${profile}"
condition_job=$(sbatch --parsable \
    --job-name=h96_C_cluster_tol \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --array="1-${task_count}%8" \
    --cpus-per-task=1 \
    --mem=32G \
    --time=03:00:00 \
    --output="${log_root}/condition_%A_%a.out" \
    --error="${log_root}/condition_%A_%a.err" \
    --export=ALL,CLUSTER_CODE_ROOT,CLUSTER_EXISTING_DATA_ROOT,CLUSTER_OUTPUT_ROOT,CLUSTER_PROFILE \
    "${code_root}/run_cluster_tolerance_task.sh")

export CLUSTER_DATA_ROOT="${data_root}"
export CLUSTER_FIGURE_ROOT="${figure_root}"
plot_job=$(sbatch --parsable \
    --job-name=h96_C_cluster_plot \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --dependency="afterok:${condition_job}" \
    --cpus-per-task=1 \
    --mem=8G \
    --time=01:00:00 \
    --output="${log_root}/plot_%j.out" \
    --error="${log_root}/plot_%j.err" \
    --export=ALL,CLUSTER_CODE_ROOT,CLUSTER_DATA_ROOT,CLUSTER_FIGURE_ROOT,CLUSTER_PROFILE \
    "${code_root}/plot_cluster_tolerance_comparison.sh")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\troot\n' > "${manifest}"
printf 'condition_extract\t%s\t%s\n' "${condition_job}" "${data_root}" >> "${manifest}"
printf 'comparison_plot\t%s\t%s\n' "${plot_job}" "${figure_root}" >> "${manifest}"
printf 'local_figure_destination\t-\t%s\n' "${local_figure_destination}" >> "${manifest}"
printf 'profile\t-\t%s\n' "${profile}" >> "${manifest}"
printf 'condition_job=%s\nplot_job=%s\nmanifest=%s\n' \
    "${condition_job}" "${plot_job}" "${manifest}"
