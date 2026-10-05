#!/bin/bash
set -euo pipefail

old_table=${1:?Usage: seed_l6_ei_fine_grid.sh OLD_TABLE NEW_RUN_ROOT SETUP_FILE}
new_run_root=${2:?Usage: seed_l6_ei_fine_grid.sh OLD_TABLE NEW_RUN_ROOT SETUP_FILE}
setup_file=${3:?Usage: seed_l6_ei_fine_grid.sh OLD_TABLE NEW_RUN_ROOT SETUP_FILE}
points_root="$new_run_root/results/points"
missing_file="$new_run_root/missing_task_ids.tsv"

mkdir -p "$points_root" "$new_run_root/results/logs"
cp "$setup_file" "$new_run_root/global_bifurcation_setup.mat"

awk -F '\t' -v OFS='\t' -v points_root="$points_root" \
    -v missing_file="$missing_file" '
  NR == 1 { next }
  {
    new_beta6_index = 2 * ($2 - 1) + 1
    new_betaEI_index = 2 * ($3 - 1) + 1
    new_task_id = (new_beta6_index - 1) * 81 + new_betaEI_index
    $1 = new_task_id
    $2 = new_beta6_index
    $3 = new_betaEI_index
    output_file = sprintf("%s/point_%05d.tsv", points_root, new_task_id)
    print > output_file
    close(output_file)
    seeded[new_task_id] = 1
  }
  END {
    for (task_id = 1; task_id <= 11421; task_id++) {
      if (!(task_id in seeded)) {
        print task_id > missing_file
      }
    }
    close(missing_file)
  }
' "$old_table"

cat > "$new_run_root/grid_config.sh" <<EOF
export L6EI_GRID_STEP=0.0025
export L6EI_BETA6_MAXIMUM=0.35
export L6EI_BETAEI_MAXIMUM=0.2
export L6EI_TOTAL_GRID_COUNT=11421
export L6EI_CHUNK_SIZE=2600
export L6EI_TASK_LIST_FILE="$missing_file"
EOF

seeded_count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' | wc -l)
missing_count=$(wc -l < "$missing_file")
if [[ "$seeded_count" -ne 2911 || "$missing_count" -ne 8510 ]]; then
  echo "Unexpected seeded/missing counts: $seeded_count seeded, $missing_count missing." >&2
  exit 1
fi
printf 'seeded_points=%d missing_points=%d total_points=%d\n' \
  "$seeded_count" "$missing_count" "$((seeded_count + missing_count))"
