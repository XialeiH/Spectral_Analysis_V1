#!/bin/bash
#SBATCH --job-name=hc10branch
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/perturbation_7_20260807/logs/hc10branch_%j.out

set -euo pipefail
module purge
module load matlab/2025b

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
CORRECTED=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle_extended_corrected_h96_ikp7_l6p5_matchedlib_20260702_010052
LONG_DATA=${RUN_ROOT}/hc10_long_horizon/hc10_long_horizon_data.mat
OUTPUT=${RUN_ROOT}/hc10_extension_branch

matlab -batch "addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); inspect_hc10_extension_branch('${CODE_ROOT}/global_bifurcation_setup.mat','${LONG_DATA}','${CORRECTED}','${OUTPUT}');"
