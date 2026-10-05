"""Collect the isolated HC-norm benchmark without querying Slurm."""
import csv
import json
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REMOTE = r'''
import csv, hashlib, json, math
from pathlib import Path

root=Path('/scratch/xh2906/librarySCI_runs/figure1b_hcnorm_residual_allfields_20260925')
expected_hashes={
    'run_figure1g2_tol5e3.m':'4f8323c173ba1c5f3d6e9463527d3c66b04ddc9f72c1002683ba4428d8d1a17e',
    'run_impl_probe.m':'c3a62a4bc0938f5c69419cfc7573d40810de4e085f3e7c24a1e3cdac6d35dc8b',
    'run_allfields.m':'8feb4cddea472195e9403df58d64d0bef4a7c765b4818cbcec268585a6afa054',
    'h96_nn_all_population_responses_symmetry.m':'96c223d8a7740b29eccf0ac07df5144be141d1d9efce32c217f57aee90c657bd',
}
errors=[]; records=[]
for name,expected in expected_hashes.items():
    path=root/'code'/name
    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
        errors.append(str(path)+': staged code hash mismatch')

for path in sorted((root/'results').glob('*.tsv')):
    try:
        with path.open(newline='') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        assert len(rows)==1
        row=rows[0]
        field=int(row['FieldHC']); repeat=int(row['Repeat'])
        is_cg=path.name.startswith('cg_dense_')
        is_dnn=path.name.startswith('paged_combined_')
        assert is_cg != is_dnn
        method='CG' if is_cg else 'DNN surrogate'
        expected_name=(f'cg_dense_{field:02d}HC_repeat_{repeat:03d}.tsv' if is_cg else
                       f'paged_combined_repeat_{repeat:03d}.tsv' if field==3 else
                       f'paged_combined_{field:02d}HC_repeat_{repeat:03d}.tsv')
        assert path.name==expected_name
        assert field in [3,4,6,8,10,20,30,40]
        expected_n=100 if field in [3,4] else 5 if is_cg and field in [30,40] else 20
        assert 1<=repeat<=expected_n
        seconds=float(row['Seconds']); iterations=int(row['Iterations'])
        residual=float(row['FinalResidual'])
        assert math.isfinite(seconds) and seconds>0
        assert math.isfinite(residual) and residual<.005
        assert 1<=iterations<=200 and str(row['Converged']).lower() in ['1','true']
        item={'FieldHC':field,'Repeat':repeat,'Method':method,'Seconds':seconds,
              'Iterations':iterations,'Converged':1,'FinalResidual':residual,
              'NNThreads':0 if is_cg else int(row['NNThreads']),
              'Variant':'cg_dense' if is_cg else str(row['Variant']),
              'Passed':1 if is_cg else int(str(row['Passed']).lower() in ['1','true'])}
        if is_dnn:
            assert row['Variant']=='paged_combined' and item['NNThreads']==8
            assert int(row['AllocatedCPUs'])==16 and item['Passed']==1
            assert iterations==int(row['ReferenceIterations'])
            assert float(row['MaxAbsolute'])<=1e-3
            assert float(row['MeanAbsolute'])<=1e-4
            assert float(row['StepMaxAbsolute'])<=1e-3
            assert float(row['StepMeanAbsolute'])<=1e-4
            assert float(row['ResidualDifference'])<=1e-6
            assert float(row['StepResidualDifference'])<=1e-6
        records.append(item)
    except Exception as exc:
        errors.append(str(path)+': '+repr(exc))

keys=[(r['FieldHC'],r['Repeat'],r['Method']) for r in records]
if len(keys)!=len(set(keys)): errors.append('Duplicate field/method/repeat')
print(json.dumps({'records':records,'errors':errors,'hashes_valid':not any('hash mismatch' in e for e in errors)}))
'''

result = subprocess.run(
    ["ssh", "-S", "/Users/xialeihuang/.ssh/torch-codex-cm", "-o", "BatchMode=yes",
     "-o", "ConnectTimeout=15", "xh2906@login.torch.hpc.nyu.edu", "python3 -"],
    input=REMOTE, text=True, capture_output=True, check=True, timeout=180)
payload = json.loads(result.stdout)
records = sorted(payload["records"], key=lambda r: (r["FieldHC"], r["Method"], r["Repeat"]))
expected = {f"{field} {method}": (100 if field in [3, 4] else
            5 if method == "CG" and field in [30, 40] else 20)
            for field in [3, 4, 6, 8, 10, 20, 30, 40]
            for method in ["CG", "DNN surrogate"]}
counts = {key: 0 for key in expected}
for row in records:
    counts[f"{row['FieldHC']} {row['Method']}"] += 1
complete = not payload["errors"] and all(counts[k] == n for k, n in expected.items())
status = {"complete": complete, "metric": "relative_HC_norm", "tolerance": 0.005,
          "expected": expected, "counts": counts, "errors": payload["errors"],
          "hashes_valid": payload["hashes_valid"]}
with (ROOT / "hcnorm_trials.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["FieldHC", "Repeat", "Method", "Seconds",
        "Iterations", "Converged", "FinalResidual", "NNThreads", "Variant", "Passed"])
    writer.writeheader()
    writer.writerows(records)
summary = []
for key, n in expected.items():
    field, method = key.split(" ", 1)
    values = [r["Seconds"] for r in records if r["FieldHC"] == int(field)
              and r["Method"] == method]
    if not values:
        continue
    sd = statistics.stdev(values) if len(values) > 1 else None
    summary.append({"FieldHC": int(field), "Method": method, "N": len(values),
                    "ExpectedN": n, "Mean": statistics.mean(values), "SD": sd,
                    "SEM": sd / len(values)**0.5 if sd is not None else None,
                    "Min": min(values), "Max": max(values)})
status["summary"] = summary
(ROOT / "status.json").write_text(json.dumps(status, indent=2) + "\n")
with (ROOT / "summary.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["FieldHC", "Method", "N", "ExpectedN",
        "Mean", "SD", "SEM", "Min", "Max"])
    writer.writeheader()
    writer.writerows(summary)
print(json.dumps({"complete": complete, "total": len(records),
                  "counts": counts, "errors": payload["errors"]}, indent=2))
if payload["errors"]:
    raise SystemExit("Validation failed; do not update the figure.")
