"""Collect current benchmark TSVs without querying Slurm or changing remote files."""
import csv
import io
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path(__file__).resolve().parent.parent
REMOTE = r'''
import csv, hashlib, json
from pathlib import Path
small=Path("/scratch/xh2906/librarySCI_runs/figure1b_tol5e3_3x3_4x4_20260921")
large=Path("/scratch/xh2906/librarySCI_runs/figure1b_tol5e3_allfields_snn15s_20260921")
expected_sha="19b51fb292c065ae59f0595aab9142bce49ab8ed342014de625a6eb490279516"
trials=[]; snn=[]; errors=[]; warnings=[]; checks=[]
for root in [small,large]:
    for path in sorted((root/"results").glob("*.tsv")):
        try:
            rows=list(csv.DictReader(path.open(),delimiter="\t"))
            assert len(rows)==1
            row=rows[0]
            for key in ["FieldHC","Repeat"]: row[key]=int(row[key])
            row["Seconds"]=float(row["Seconds"])
            assert row["Seconds"]>0
            if row["Method"]=="SNN":
                for key in ["SimulationMs","dtMs","Updates","CPUs"]: row[key]=float(row[key])
                assert (row["FieldHC"],row["SimulationMs"],row["dtMs"],row["Updates"],row["CPUs"])==(3,15000,0.1,150000,16)
                snn.append(row)
            else:
                row["Iterations"]=int(row["Iterations"])
                row["FinalResidual"]=float(row["FinalResidual"])
                row["Converged"]=int(row["Converged"].lower() in ("1","true"))
                assert row["Converged"]==1 and row["FinalResidual"]<.005 and row["Iterations"]<=200
                row["Checksum"]=float(row["Checksum"])
                trials.append(row)
        except Exception as exc:
            errors.append(str(path)+": "+repr(exc))
    for path in sorted((root/"logs").glob("*.err")):
        if path.stat().st_size:
            warnings.append({"path":str(path),"tail":path.read_text(errors="replace")[-3000:]})
for i in range(1,6):
    trial=large/("snn_trial_%03d"%i)
    driver=trial/"Paper2_Fig7Comp_NW_LDE.m"
    if not driver.exists(): continue
    sha=hashlib.sha256(driver.read_bytes()).hexdigest()
    if sha!=expected_sha: errors.append(str(driver)+": SHA256 mismatch")
    folder=trial/"Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData"
    segments=list(folder.glob("DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat"))
    final=folder/"NWSimulationPix_0.0deg_NewSmear.mat"
    valid=(sha==expected_sha and len(segments)==15 and all(p.stat().st_size>0 for p in segments) and final.exists() and final.stat().st_size>0)
    checks.append({"repeat":i,"segments":len(segments),"driver_sha":sha,"valid":valid})
    if any(r["Repeat"]==i for r in snn) and not valid: errors.append(str(trial)+": completed SNN failed artifact validation")
keys=[(r["FieldHC"],r["Method"],r["Repeat"]) for r in trials+snn]
if len(keys)!=len(set(keys)): errors.append("Duplicate field/method/repeat")
print(json.dumps({"trials":trials,"snn":snn,"validation_errors":errors,"log_warnings":warnings,"snn_checks":checks,"snn_validated":len(snn)==5 and len(checks)==5 and all(c["valid"] for c in checks)}))
'''
result = subprocess.run(
    ["ssh", "-S", "/Users/xialeihuang/.ssh/torch-codex-cm",
     "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
     "xh2906@login.torch.hpc.nyu.edu", "python3 -"],
    input=REMOTE, text=True, capture_output=True, check=True, timeout=180)
payload = json.loads(result.stdout)
trial_fields = ["FieldHC","Repeat","Method","Iterations","Seconds","Converged","FinalResidual","Checksum"]
snn_fields = ["FieldHC","Repeat","Method","Seconds","SimulationMs","dtMs","Updates","CPUs","Protocol"]
for name, rows, fields in [("current_trials.csv",payload["trials"],trial_fields),
                            ("current_snn.csv",payload["snn"],snn_fields)]:
    with (ROOT/name).open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
status={k:v for k,v in payload.items() if k not in ("trials","snn")}
status["counts"]={}
summary=[]
for field in [3,4,6,8,10,20,30,40]:
    for method in ["CG","DNN surrogate","SNN"]:
        rows=[r for r in payload["trials"]+payload["snn"] if r["FieldHC"]==field and r["Method"]==method]
        if not rows: continue
        values=[r["Seconds"] for r in rows]
        sd=statistics.stdev(values) if len(values)>1 else None
        row={"FieldHC":field,"Method":method,"N":len(values),
             "Mean":statistics.mean(values),"SD":sd,
             "SEM":sd/len(values)**.5 if sd is not None else None,
             "Min":min(values),"Max":max(values)}
        summary.append(row)
        status["counts"][f"{field}x{field} {method}"]=len(values)
status["summary"]=summary
(ROOT/"current_status.json").write_text(json.dumps(status,indent=2)+"\n")
with (ROOT/"current_summary.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=["FieldHC","Method","N","Mean","SD","SEM","Min","Max"])
    writer.writeheader(); writer.writerows(summary)
print(json.dumps(status,indent=2))
if status["validation_errors"]:
    raise SystemExit("Validation failed; do not update the figure.")
