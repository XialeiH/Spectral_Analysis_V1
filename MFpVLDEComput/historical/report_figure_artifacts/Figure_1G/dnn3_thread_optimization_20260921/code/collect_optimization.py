"""Fetch isolated DNN optimization results; never substitute target-derived times."""
import csv
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path(__file__).resolve().parent.parent
REMOTE = r'''
import csv,json
from pathlib import Path
r=Path("/scratch/xh2906/librarySCI_runs/figure1b_dnn3_thread_optimization_20260921")
def rows(pattern):
    out=[]
    for p in sorted(r.glob(pattern)):
        with p.open() as f: out.extend(csv.DictReader(f,delimiter="\t"))
    return out
selected=r/"selected_threads.txt"
summary=r/"probe_summary.csv"
warnings=[{"path":str(p),"text":p.read_text(errors="replace")[-3000:]} for p in (r/"logs").glob("*.err") if p.stat().st_size]
print(json.dumps({"probe":rows("probe/*.tsv"),"validation":rows("validation/dnn_*.tsv"),
    "checks":rows("validation/check_*.tsv"),
    "selected_threads":int(selected.read_text()) if selected.exists() else None,
    "probe_summary":summary.read_text() if summary.exists() else None,
    "log_warnings":warnings}))
'''
response=subprocess.run(["ssh","-S","/Users/xialeihuang/.ssh/torch-codex-cm",
    "-o","BatchMode=yes","-o","ConnectTimeout=15","xh2906@login.torch.hpc.nyu.edu","python3 -"],
    input=REMOTE,text=True,capture_output=True,check=True,timeout=180)
data=json.loads(response.stdout)
trial_fields=["FieldHC","Repeat","Method","Iterations","Seconds","Converged",
    "FinalResidual","Checksum","NNThreads","AllocatedCPUs","Node"]
check_fields=["Repeat","NNThreads","MaxAbs","MeanAbs","SameIterations","ResidualDifference","Passed"]
for name,rows,fields in [("probe_trials.csv",data["probe"],trial_fields),
    ("optimized_trials.csv",data["validation"],trial_fields),
    ("equivalence_checks.csv",data["checks"],check_fields)]:
    with (ROOT/name).open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
if data["probe_summary"] is not None:
    (ROOT/"probe_summary.csv").write_text(data["probe_summary"])
status={"probe_count":len(data["probe"]),"validation_count":len(data["validation"]),
    "check_count":len(data["checks"]),"selected_threads":data["selected_threads"],
    "log_warnings":data["log_warnings"],"complete":False,"ready_to_plot":False}
yes=lambda value: str(value).lower() in ("1","true")
if len(data["validation"])==100 and len(data["checks"])==100:
    trials=sorted(data["validation"],key=lambda r:int(r["Repeat"]))
    checks=sorted(data["checks"],key=lambda r:int(r["Repeat"]))
    assert [int(r["Repeat"]) for r in trials]==list(range(1,101))
    assert [int(r["Repeat"]) for r in checks]==list(range(1,101))
    assert all(yes(r["Converged"]) and float(r["FinalResidual"])<.005 and int(r["FieldHC"])==3
               and int(r["NNThreads"])==data["selected_threads"] and int(r["AllocatedCPUs"])==16 for r in trials)
    assert all(yes(r["Passed"]) and yes(r["SameIterations"]) and float(r["MaxAbs"])<=.001
               and float(r["MeanAbs"])<=.0001 and float(r["ResidualDifference"])<=.000001 for r in checks)
    times=[float(r["Seconds"]) for r in trials]
    assert all(t>0 for t in times)
    mean=statistics.mean(times); sd=statistics.stdev(times)
    status.update(complete=True,mean=mean,sd=sd,sem=sd/10,min=min(times),max=max(times),
        cg_dnn=27.36131558/mean,snn_dnn=979.9515794754/mean,
        target_100x_met=mean<=27.36131558/100,
        ready_to_plot=mean<.53661261 and mean<.54671737,
        max_state_difference=max(float(r["MaxAbs"]) for r in checks))
(ROOT/"optimization_status.json").write_text(json.dumps(status,indent=2)+"\n")
print(json.dumps(status,indent=2))
