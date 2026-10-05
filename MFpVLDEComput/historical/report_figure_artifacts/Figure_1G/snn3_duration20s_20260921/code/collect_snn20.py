"""Collect measured20s SNN full-process times without changing figure inputs."""
import csv
import json
from pathlib import Path
import statistics
import subprocess

ROOT=Path(__file__).resolve().parent.parent
REMOTE=r'''
import csv,hashlib,json
from pathlib import Path
r=Path("/scratch/xh2906/librarySCI_runs/figure1b_snn3_duration20s_20260921")
sha="b3958be55c41dba44f674de73e3dfb8bf37e2b552c7f4b145b824872197f4c6b"
rows=[]; validation=[]; errors=[]
for p in sorted((r/"results").glob("snn_paper2_20s_03HC_repeat_*.tsv")):
    data=list(csv.DictReader(p.open(),delimiter="\t"))
    assert len(data)==1
    d=data[0]; i=int(d["Repeat"])
    assert (int(d["FieldHC"]),float(d["SimulationMs"]),float(d["dtMs"]),int(d["Updates"]),int(d["CPUs"]))==(3,20000,.1,200000,16)
    assert float(d["Seconds"])>0
    trial=r/("snn_trial_%03d"%i)
    driver=trial/"Paper2_Fig7Comp_NW_LDE.m"
    folder=trial/"Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData"
    found=list(folder.glob("DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat"))
    segments=[folder/("DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_%ds_NewSmear.mat"%j) for j in range(1,21)]
    final=folder/"NWSimulationPix_0.0deg_NewSmear.mat"
    digest=hashlib.sha256(driver.read_bytes()).hexdigest()
    valid=digest==sha and len(found)==20 and all(f.exists() and f.stat().st_size>0 for f in segments+[final])
    validation.append({"repeat":i,"segments":len(found),"sha256":digest,"valid":valid})
    if not valid: errors.append(str(trial)+": invalid completed artifacts")
    rows.append(d)
warnings=[{"path":str(p),"tail":p.read_text(errors="replace")[-3000:]} for p in (r/"logs").glob("*.err") if p.stat().st_size]
print(json.dumps({"trials":rows,"validation":validation,"errors":errors,"log_warnings":warnings}))
'''
r=subprocess.run(["ssh","-S","/Users/xialeihuang/.ssh/torch-codex-cm",
    "-o","BatchMode=yes","-o","ConnectTimeout=15","xh2906@login.torch.hpc.nyu.edu","python3 -"],
    input=REMOTE,text=True,capture_output=True,check=True,timeout=180)
data=json.loads(r.stdout)
fields=["FieldHC","Repeat","Method","Seconds","SimulationMs","dtMs","Updates","CPUs","Protocol"]
with (ROOT/"trials.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(data["trials"])
times=[float(d["Seconds"]) for d in data["trials"]]
data["complete"]=len(times)==5 and not data["errors"] and all(c["valid"] for c in data["validation"])
if data["complete"]:
    assert sorted(int(d["Repeat"]) for d in data["trials"])==list(range(1,6))
    mean=statistics.mean(times); sd=statistics.stdev(times)
    data["summary"]={"N":5,"Mean":mean,"SD":sd,"SEM":sd/5**.5,
        "Min":min(times),"Max":max(times),"RatioTo15sMean":mean/979.9515794754}
    with (ROOT/"summary.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data["summary"])); w.writeheader(); w.writerow(data["summary"])
(ROOT/"status.json").write_text(json.dumps(data,indent=2)+"\n")
print(json.dumps(data,indent=2))
if data["errors"]: raise SystemExit("SNN20 validation failed.")
