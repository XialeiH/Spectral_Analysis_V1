#!/usr/bin/env bash
#SBATCH --job-name=l6ns_atlas_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/atlas_aggregate_%j.out
#SBATCH --error=logs/atlas_aggregate_%j.err

set -euo pipefail
module load matlab/2025b
cd "${L6NS_ATLAS_ROOT}"
matlab -nodisplay -nosplash -batch "addpath('${L6NS_ATLAS_ROOT}','-begin'); aggregate_l6ns_branch_atlas_phase1('${L6NS_ATLAS_ROOT}')"
