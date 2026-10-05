#!/usr/bin/env bash
set -euo pipefail

module purge
module load matlab/2025b

matlab -batch "addpath('${CLUSTER_CODE_ROOT}'); run_jfull_c_cluster_tolerance_task();"
