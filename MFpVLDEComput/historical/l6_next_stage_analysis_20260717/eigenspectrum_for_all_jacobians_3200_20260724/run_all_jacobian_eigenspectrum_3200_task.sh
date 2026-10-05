#!/bin/bash
set -euo pipefail

module load matlab/2025b
cd "$ALLJ_CODE_ROOT"
matlab -nodisplay -nosplash -batch "run_all_jacobian_eigenspectrum_3200_task"
