#!/bin/bash
#SBATCH --job-name=hc5ftall
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-5%2
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/hc5ftall_%A_%a.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
SOURCE=${RUN_ROOT}/outputs/7.6_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC5_data.mat
OUTPUT_ROOT=${RUN_ROOT}/hc5_finite_time_all_snapshots

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); compare_hc5_snapshot_finite_time_task(str2double(getenv('SLURM_ARRAY_TASK_ID')),'${CODE_ROOT}/global_bifurcation_setup.mat','${SOURCE}','${OUTPUT_ROOT}');"
