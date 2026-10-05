#!/usr/bin/env bash
set -euo pipefail

module purge
module load matlab/2025b

matlab -batch "addpath('${FIXED_CLUSTER_CODE_ROOT}'); run_jfull_c_fixed_mode_count_task();"
