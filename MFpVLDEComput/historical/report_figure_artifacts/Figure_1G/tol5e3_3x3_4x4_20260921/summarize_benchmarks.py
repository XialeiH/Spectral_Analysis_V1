"""Summarize measured iteration-only timings; never infer missing trials."""
import argparse
import csv
import io
import json
from pathlib import Path
import statistics
import tarfile

ROOT = Path(__file__).resolve().parent

def read_trials(directory):
    return [row for path in sorted(directory.glob("*.tsv"))
            for row in csv.DictReader(path.open(), delimiter="\t")]

def summarize(rows, repeats):
    records = []
    assert len(rows) == 4 * repeats
    for field in (3, 4):
        for method in ("CG", "DNN surrogate"):
            group = [r for r in rows if int(r["FieldHC"]) == field and r["Method"] == method]
            assert sorted(int(r["Repeat"]) for r in group) == list(range(1, repeats + 1))
            assert all(r["Converged"] == "1" and float(r["FinalResidual"]) < 5e-3 for r in group)
            seconds = [float(r["Seconds"]) for r in group]
            iterations = [int(r["Iterations"]) for r in group]
            residuals = [float(r["FinalResidual"]) for r in group]
            records.append(dict(field_HC=field, method=method, n=len(group),
                mean_seconds=statistics.mean(seconds), sd_seconds=statistics.stdev(seconds),
                mean_iterations=statistics.mean(iterations), min_iterations=min(iterations),
                max_iterations=max(iterations), min_residual=min(residuals), max_residual=max(residuals)))
    return records

def comparison(records):
    out = []
    for field in (3, 4):
        cg = next(r for r in records if r["field_HC"] == field and r["method"] == "CG")
        dnn = next(r for r in records if r["field_HC"] == field and r["method"] == "DNN surrogate")
        out.append(dict(field_HC=field, mean_CG_seconds=cg["mean_seconds"],
            mean_DNN_seconds=dnn["mean_seconds"],
            CG_over_DNN=cg["mean_seconds"]/dnn["mean_seconds"]))
    return out

parser = argparse.ArgumentParser()
parser.add_argument("--local-only", action="store_true")
args = parser.parse_args()
local = summarize(read_trials(ROOT/"local_smoke/output/results"), 3)
report = dict(relative_step_tolerance=5e-3, consecutive_passes=3,
    maximum_iterations=200, relaxation=0.33, condition=dict(angle_deg=0, contrast=100),
    speedup_definition="ratio of arithmetic mean times",
    local_hardware="Apple M3 Pro, 12 physical/logical CPUs, 36 GiB RAM",
    local_thread_limit=12, local_trials=local, local_speedups=comparison(local),
    timing_scope="CG/DNN: iteration loop including convergence checks; excludes load/preparation. SNN unchanged.")
if not args.local_only:
    torch = summarize(read_trials(ROOT/"torch_snapshot/results"), 100)
    report.update(torch_trials=torch, torch_speedups=comparison(torch),
        torch_CPUs_per_trial=16, torch_GiB_per_trial=32)
    with tarfile.open(ROOT/"original_strict_trials.tar.gz") as archive:
        original = [r for item in archive.getmembers() if item.isfile()
                    for r in csv.DictReader(io.TextIOWrapper(archive.extractfile(item)), delimiter="\t")]
    snn = [float(r["Seconds"]) for r in original if r["Method"] == "SNN" and int(r["FieldHC"]) == 3]
    assert len(snn) == 5
    report["unchanged_SNN_3x3_mean_seconds"] = statistics.mean(snn)
    report["SNN_over_CG_3x3"] = statistics.mean(snn)/report["torch_speedups"][0]["mean_CG_seconds"]
name = "local_smoke_summary.json" if args.local_only else "benchmark_summary.json"
(ROOT/name).write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps(report, indent=2))

