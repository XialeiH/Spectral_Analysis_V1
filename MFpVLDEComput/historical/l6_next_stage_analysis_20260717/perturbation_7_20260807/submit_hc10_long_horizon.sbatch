#!/bin/bash
#SBATCH --job-name=hc10long
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/hc10long_%j.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
SOURCE=${RUN_ROOT}/outputs/7.7_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC10_data.mat
OUTPUT=${RUN_ROOT}/hc10_long_horizon

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); run_hc10_long_horizon('${CODE_ROOT}/global_bifurcation_setup.mat','${SOURCE}','${OUTPUT}');"
