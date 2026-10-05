#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 REMOTE_CODE_ROOT REMOTE_OUTPUT_ROOT" >&2
  exit 2
fi

CODE_ROOT=$1
OUTPUT_ROOT=$2
BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
L6_SOURCE=$BASE/l6_mechanism_derivative_clamp_20260628_033413
INHIBITION_SOURCE=$BASE/l6_celltype_controls_20260628_034009
L6_ATLAS=/scratch/xh2906/librarySCI_runs/h96_l6_multiangle_top4_atlas_20260720_115904
INHIBITION_ATLAS=/scratch/xh2906/librarySCI_runs/h96_inhibition_multiangle_top4_atlas_20260720_125113

mkdir -p "$OUTPUT_ROOT" "$OUTPUT_ROOT/logs"
rm -rf "$OUTPUT_ROOT/L6" "$OUTPUT_ROOT/Inhibition"
rm -f "$OUTPUT_ROOT"/manifest_*.tsv "$OUTPUT_ROOT"/boundary_*.tsv "$OUTPUT_ROOT"/logs/*

export ALLJ_CODE_ROOT="$CODE_ROOT"
export ALLJ_OUTPUT_ROOT="$OUTPUT_ROOT"
export ALLJ_L6_SOURCE_ROOT="$L6_SOURCE"
export ALLJ_INHIBITION_SOURCE_ROOT="$INHIBITION_SOURCE"
export ALLJ_L6_ATLAS_ROOT="$L6_ATLAS"
export ALLJ_INHIBITION_ATLAS_ROOT="$INHIBITION_ATLAS"

job_id=$(sbatch --parsable \
  --job-name=allJ3200 \
  --account=torch_pr_482_general \
  --partition=cpu_short \
  --qos=cpu48 \
  --array=1-8%8 \
  --cpus-per-task=4 \
  --mem=32G \
  --time=02:00:00 \
  --output="$OUTPUT_ROOT/logs/allJ3200_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/allJ3200_%A_%a.err" \
  --export=ALL \
  "$CODE_ROOT/run_all_jacobian_eigenspectrum_3200_task.sh")

manifest="$OUTPUT_ROOT/submission_manifest.tsv"
printf 'job_id\tcode_root\toutput_root\n%s\t%s\t%s\n' \
  "$job_id" "$CODE_ROOT" "$OUTPUT_ROOT" > "$manifest"
printf 'job_id=%s\nmanifest=%s\n' "$job_id" "$manifest"
