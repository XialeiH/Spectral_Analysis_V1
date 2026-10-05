#!/bin/bash
set -euo pipefail

module load matlab/2025b
cd "$JI_CODE_ROOT"
matlab -nodisplay -nosplash -batch "run_ji_component_structure_task"
