#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 RUN_ROOT EXISTING_DATA_ROOT [PROFILE]" >&2
    exit 2
fi

run_root=$1
existing_data_root=$2
profile=${3:-initial}
case "${profile}" in
    initial|revised) ;;
    *) echo "Unsupported profile: ${profile}" >&2; exit 2 ;;
esac
code_root="${run_root}/code"
data_root="${run_root}/data"
figure_root="${run_root}/figures"
log_root="${run_root}/logs"
mkdir -p "${data_root}" "${figure_root}" "${log_root}"

export FIXED_CLUSTER_CODE_ROOT="${code_root}"
export FIXED_CLUSTER_EXISTING_DATA_ROOT="${existing_data_root}"
export FIXED_CLUSTER_OUTPUT_ROOT="${data_root}"
export FIXED_CLUSTER_PROFILE="${profile}"
condition_job=$(sbatch --parsable \
    --job-name=h96_C_fixed_cluster \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --array=1-8%8 \
    --cpus-per-task=1 \
    --mem=8G \
    --time=01:00:00 \
    --output="${log_root}/condition_%A_%a.out" \
    --error="${log_root}/condition_%A_%a.err" \
    --export=ALL,FIXED_CLUSTER_CODE_ROOT,FIXED_CLUSTER_EXISTING_DATA_ROOT,FIXED_CLUSTER_OUTPUT_ROOT,FIXED_CLUSTER_PROFILE \
    "${code_root}/run_fixed_mode_count_task.sh")

export FIXED_CLUSTER_DATA_ROOT="${data_root}"
export FIXED_CLUSTER_FIGURE_ROOT="${figure_root}"
plot_job=$(sbatch --parsable \
    --job-name=h96_C_fixed_plot \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --dependency="afterok:${condition_job}" \
    --cpus-per-task=1 \
    --mem=8G \
    --time=01:00:00 \
    --output="${log_root}/plot_%j.out" \
    --error="${log_root}/plot_%j.err" \
    --export=ALL,FIXED_CLUSTER_CODE_ROOT,FIXED_CLUSTER_DATA_ROOT,FIXED_CLUSTER_FIGURE_ROOT,FIXED_CLUSTER_PROFILE \
    "${code_root}/plot_fixed_mode_count_comparison.sh")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\troot\n' > "${manifest}"
printf 'condition_extract\t%s\t%s\n' "${condition_job}" "${data_root}" >> "${manifest}"
printf 'comparison_plot\t%s\t%s\n' "${plot_job}" "${figure_root}" >> "${manifest}"
printf 'profile\t-\t%s\n' "${profile}" >> "${manifest}"
printf 'condition_job=%s\nplot_job=%s\nmanifest=%s\n' \
    "${condition_job}" "${plot_job}" "${manifest}"
