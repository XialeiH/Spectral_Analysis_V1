#!/usr/bin/env bash
set -euo pipefail

module purge
module load matlab/2025b

matlab -batch "addpath('${FULLMAP_CODE_ROOT}'); run_full_eigenspectrum_maps_condition_task();"
