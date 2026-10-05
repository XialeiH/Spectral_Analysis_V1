#!/usr/bin/env bash
#SBATCH --job-name=l6ns_atlas
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=03:30:00
#SBATCH --output=logs/atlas_%A_%a.out
#SBATCH --error=logs/atlas_%A_%a.err

set -euo pipefail
module load matlab/2025b
cd "${L6NS_ATLAS_ROOT}"
matlab -nodisplay -nosplash -batch "addpath('${L6NS_ATLAS_ROOT}','-begin'); addpath('${L6NS_ATLAS_ROOT}/Utils','-begin'); addpath('${L6NS_ATLAS_ROOT}/real_tuning_l6_i_20260722/mechanism_analysis_20260722/runtime_h96','-begin'); run_l6ns_branch_atlas_task"
