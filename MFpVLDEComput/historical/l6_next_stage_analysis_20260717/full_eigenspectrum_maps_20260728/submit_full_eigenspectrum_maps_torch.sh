#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 RUN_ROOT SOURCE_ROOT EIGEN_CACHE_ROOT FIGURE_ROOT" >&2
    exit 2
fi

run_root=$1
source_root=$2
eigen_cache_root=$3
figure_root=$4
code_root="${run_root}/code"
data_root="${run_root}/results"
log_root="${run_root}/logs"
mkdir -p "${data_root}" "${figure_root}" "${log_root}"

export FULLMAP_CODE_ROOT="${code_root}"
export FULLMAP_SOURCE_ROOT="${source_root}"
export FULLMAP_EIGEN_CACHE_ROOT="${eigen_cache_root}"
export FULLMAP_OUTPUT_ROOT="${data_root}"

condition_job=$(sbatch --parsable \
    --job-name=h96_fullmap_modes \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --array=1-16%8 \
    --cpus-per-task=1 \
    --mem=32G \
    --time=03:00:00 \
    --output="${log_root}/condition_%A_%a.out" \
    --error="${log_root}/condition_%A_%a.err" \
    --export=ALL,FULLMAP_CODE_ROOT,FULLMAP_SOURCE_ROOT,FULLMAP_EIGEN_CACHE_ROOT,FULLMAP_OUTPUT_ROOT \
    "${code_root}/run_full_eigenspectrum_maps_condition_task.sh")

export FULLMAP_DATA_ROOT="${data_root}/data"
export FULLMAP_FIGURE_ROOT="${figure_root}"
plot_job=$(sbatch --parsable \
    --job-name=h96_fullmap_plot \
    --account=torch_pr_482_general \
    --partition=cpu_short \
    --qos=cpu48 \
    --dependency="afterok:${condition_job}" \
    --cpus-per-task=1 \
    --mem=32G \
    --time=03:00:00 \
    --output="${log_root}/plot_%j.out" \
    --error="${log_root}/plot_%j.err" \
    --export=ALL,FULLMAP_CODE_ROOT,FULLMAP_DATA_ROOT,FULLMAP_FIGURE_ROOT \
    "${code_root}/plot_full_eigenspectrum_map_atlases.sh")

manifest="${run_root}/submission_manifest.tsv"
printf 'stage\tjob_id\troot\n' > "${manifest}"
printf 'condition_extract\t%s\t%s\n' "${condition_job}" "${data_root}" >> "${manifest}"
printf 'atlas_plot\t%s\t%s\n' "${plot_job}" "${figure_root}" >> "${manifest}"
printf 'condition_job=%s\nplot_job=%s\nmanifest=%s\n' \
    "${condition_job}" "${plot_job}" "${manifest}"
