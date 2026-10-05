#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 RUN_ROOT EXISTING_DATA_ROOT" >&2
    exit 2
fi

run_root=$1
existing_data_root=$2
code_root="${run_root}/code"
data_root="${run_root}/data"
figure_root="${run_root}/figures"
log_root="${run_root}/logs"
mkdir -p "${data_root}" "${figure_root}" "${log_root}"

export REVISED_CLUSTER_CODE_ROOT="${code_root}"
export REVISED_CLUSTER_EXISTING_DATA_ROOT="${existing_data_root}"
export REVISED_CLUSTER_OUTPUT_ROOT="${data_root}"
condition_job=$(sbatch --parsable \
    --job-name=h96_Jfull_cluster \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --array=1-16%8 \
    --cpus-per-task=1 \
    --mem=8G \
    --time=01:00:00 \
    --output="${log_root}/condition_%A_%a.out" \
    --error="${log_root}/condition_%A_%a.err" \
    --export=ALL,REVISED_CLUSTER_CODE_ROOT,REVISED_CLUSTER_EXISTING_DATA_ROOT,REVISED_CLUSTER_OUTPUT_ROOT \
    "${code_root}/run_revised_jfull_cluster_task.sh")

export FULLMAP_CODE_ROOT="${code_root}"
export FULLMAP_DATA_ROOT="${data_root}"
export FULLMAP_FIGURE_ROOT="${figure_root}"
export FULLMAP_JACOBIAN_FILTER="J_full"
plot_job=$(sbatch --parsable \
    --job-name=h96_Jfull_atlases \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --dependency="afterok:${condition_job}" \
    --cpus-per-task=1 \
    --mem=8G \
    --time=03:00:00 \
    --output="${log_root}/plot_%j.out" \
    --error="${log_root}/plot_%j.err" \
    --export=ALL,FULLMAP_CODE_ROOT,FULLMAP_DATA_ROOT,FULLMAP_FIGURE_ROOT,FULLMAP_JACOBIAN_FILTER \
    "${code_root}/plot_full_eigenspectrum_map_atlases.sh")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\troot\n' > "${manifest}"
printf 'cluster_revision\t%s\t%s\n' "${condition_job}" "${data_root}" >> "${manifest}"
printf 'J_full_atlases\t%s\t%s\n' "${plot_job}" "${figure_root}" >> "${manifest}"
printf 'condition_job=%s\nplot_job=%s\nmanifest=%s\n' \
    "${condition_job}" "${plot_job}" "${manifest}"
