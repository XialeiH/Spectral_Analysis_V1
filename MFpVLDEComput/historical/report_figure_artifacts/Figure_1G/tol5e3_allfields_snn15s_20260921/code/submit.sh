#!/bin/bash
set -euo pipefail
export OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1b_tol5e3_allfields_snn15s_20260921
export SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_field_scaling_20260825_001500
export BASE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829
export EXACT_ROOT="$BASE_ROOT/snn_paper2_original_exact_20260901_160003"
export DRIVER_SHA=19b51fb292c065ae59f0595aab9142bce49ab8ed342014de625a6eb490279516
test ! -e "$OUTPUT_ROOT/submission.txt"
mkdir -p "$OUTPUT_ROOT/results" "$OUTPUT_ROOT/logs"
echo "11981ea1ae5e8927e171930770ff7858997b4c78d01755fc0f294b7f21fccd2b  $BASE_ROOT/snn_repo/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
echo "$DRIVER_SHA  $OUTPUT_ROOT/code/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
python3 -c 'import pathlib,os; b=pathlib.Path(os.environ["BASE_ROOT"])/"snn_repo/Paper2_Fig7Comp_NW_LDE.m"; p=pathlib.Path(os.environ["OUTPUT_ROOT"])/"code/Paper2_Fig7Comp_NW_LDE.m"; s=b.read_bytes(); old=b"T = 10000; dt = 0.1; TPar = 10;"; assert s.count(old)==1; assert p.read_bytes().rstrip(b"\n")==s.replace(old,b"T = 15000; dt = 0.1; TPar = 15;").rstrip(b"\n"); print("Verified: only duration/partition-count line and final newline changed.")'
cp "$BASE_ROOT/code/h96_nn_all_population_responses_symmetry.m" "$OUTPUT_ROOT/code/"
cp "$BASE_ROOT/code/h96_nn_all_population_responses_symmetry_nosparse.m" "$OUTPUT_ROOT/code/"
for field in 06 08 10 20 30 40; do
    test -s "$SOURCE_ROOT/bundles/field_${field}HC.mat"
done
bash -n "$OUTPUT_ROOT/code/large_fields.sbatch"
bash -n "$OUTPUT_ROOT/code/snn_15s.sbatch"
sbatch --parsable --export=ALL --output="$OUTPUT_ROOT/logs/snn_%A_%a.out" --error="$OUTPUT_ROOT/logs/snn_%A_%a.err" "$OUTPUT_ROOT/code/snn_15s.sbatch" | tee "$OUTPUT_ROOT/submission.txt"
sbatch --parsable --export=ALL --output="$OUTPUT_ROOT/logs/large_%A_%a.out" --error="$OUTPUT_ROOT/logs/large_%A_%a.err" "$OUTPUT_ROOT/code/large_fields.sbatch" | tee -a "$OUTPUT_ROOT/submission.txt"
