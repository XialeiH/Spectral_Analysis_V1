#!/bin/bash
#SBATCH --job-name=hcboundfig
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-6%6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/hcboundfig_%A_%a.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
OUTPUT_ROOT=${RUN_ROOT}/fixed_direction_boundary_figures

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); run_fixed_direction_perturbation_figure_task(str2double(getenv('SLURM_ARRAY_TASK_ID')),'${CODE_ROOT}/global_bifurcation_setup.mat','${OUTPUT_ROOT}');"
