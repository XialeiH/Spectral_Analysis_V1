#!/bin/bash
set -euo pipefail

module purge
module load matlab/2025b

cd "${NZ_CODE_ROOT}"
matlab -nodisplay -nosplash -batch "run_near_zero_condition_task([], '${NZ_SOURCE_ROOT}', '${NZ_OUTPUT_ROOT}');"
