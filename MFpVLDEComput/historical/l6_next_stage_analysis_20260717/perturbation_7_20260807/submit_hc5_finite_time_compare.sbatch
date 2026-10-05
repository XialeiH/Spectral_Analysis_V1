#!/bin/bash
#SBATCH --job-name=hc5ftcompare
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/hc5ftcompare_%j.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
OUTPUT_ROOT=${RUN_ROOT}/hc5_finite_time_comparison

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); compare_hc5_snapshot_finite_time_transient('${CODE_ROOT}/global_bifurcation_setup.mat','${OUTPUT_ROOT}');"
