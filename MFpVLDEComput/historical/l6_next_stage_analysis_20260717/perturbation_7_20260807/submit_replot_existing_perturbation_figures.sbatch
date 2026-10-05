#!/bin/bash
#SBATCH --job-name=replotold7
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/replotold7_%j.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
DATA_ROOT=${RUN_ROOT}/outputs
OUTPUT_ROOT=${RUN_ROOT}/renumbered_existing_figures

matlab -batch "addpath('${CODE_ROOT}','-begin'); replot_existing_perturbation_figures('${DATA_ROOT}','${OUTPUT_ROOT}');"
