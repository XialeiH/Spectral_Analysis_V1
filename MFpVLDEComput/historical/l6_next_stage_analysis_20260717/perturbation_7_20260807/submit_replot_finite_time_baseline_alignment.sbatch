#!/bin/bash
#SBATCH --job-name=fig7ftalign
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/fig7ftalign_%j.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
FINITE_ROOT=${RUN_ROOT}/all_figures_finite_time_rows
OUTPUT_ROOT=${RUN_ROOT}/all_figures_finite_time_baseline_alignment_continuous_pdfs

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); replot_all_perturbation_figures_with_finite_time('${RUN_ROOT}','${FINITE_ROOT}','${OUTPUT_ROOT}');"
