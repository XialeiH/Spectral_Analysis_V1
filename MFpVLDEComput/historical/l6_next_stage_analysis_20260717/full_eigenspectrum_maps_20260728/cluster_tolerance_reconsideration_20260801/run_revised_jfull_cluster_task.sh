#!/usr/bin/env bash
set -euo pipefail

module purge
module load matlab/2025b

matlab -batch "addpath('${REVISED_CLUSTER_CODE_ROOT}'); run_revised_jfull_cluster_condition_task();"
