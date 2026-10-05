"""Fetch actual small result files; no Slurm polling and no figure changes."""
import csv
import io
import json
from pathlib import Path
import statistics
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
REMOTE='/scratch/xh2906/librarySCI_runs/figure1b_snn3_thread_optimization_20260921'


def summarize(rows):
    if not rows: return {'n':0}
    t=[float(r['Seconds']) for r in rows]
    sd=statistics.stdev(t) if len(t)>1 else None
    return dict(n=len(t),mean=statistics.mean(t),sd=sd,
                sem=sd/len(t)**.5 if sd is not None else None,min=min(t),max=max(t))


def analyze(root):
    rows=[]
    for p in sorted((root/'results').glob('*.tsv')):
        with p.open() as stream:
            values=list(csv.DictReader(stream,delimiter='\t'))
        assert len(values)==1
        r=values[0]
        assert (int(r['FieldHC']),int(r['SimulationMs']),int(r['Updates']),int(r['Segments']),int(r['CPUs']))==(3,20000,200000,20,16)
        assert r['DriverSHA']=='b3958be55c41dba44f674de73e3dfb8bf37e2b552c7f4b145b824872197f4c6b'
        assert float(r['Seconds'])>0 and float(r['dtMs'])==.1
        rows.extend(values)
    assert len({(r['Phase'],r['Threads'],r['Repeat']) for r in rows})==len(rows)
    probe={str(t):summarize([r for r in rows if r['Phase']=='probe' and int(r['Threads'])==t]) for t in [1,2,4,8,16]}
    path=root/'selection.json'
    selection=json.loads(path.read_text()) if path.exists() else None
    pairs=[json.loads(p.read_text()) for p in sorted((root/'results').glob('pair_*.json'))]
    failures=[str(p) for p in (root/'results').glob('*.failure.json')]
    candidate=[r for r in rows if r['Phase']=='candidate']
    control=[r for r in rows if r['Phase']=='control']
    a,b=summarize(candidate),summarize(control)
    complete=False; ready=False
    if selection and selection['threads']==16:
        complete=all(v['n']==3 for v in probe.values()) and not failures
    elif selection and len(pairs)==5 and len(candidate)==len(control)==5:
        assert sorted(int(r['Repeat']) for r in candidate)==list(range(1,6))
        assert sorted(int(r['Repeat']) for r in control)==list(range(1,6))
        for c in candidate:
            ref=next(r for r in control if r['Repeat']==c['Repeat'])
            assert c['Node']==ref['Node'] and int(ref['Threads'])==16
            assert int(c['Threads'])==selection['threads']
        complete=True
        ready=not failures and all(p['passed'] for p in pairs) and a['mean']<b['mean']
    report=dict(probe=probe,selection=selection,candidate=a,control=b,pair_checks=pairs,
                failures=failures,complete=complete,ready_to_publish=ready)
    if candidate and control: report['control_over_candidate']=b['mean']/a['mean']
    if candidate: report['historical20s_over_candidate']=1436.5994801044/a['mean']
    (root/'status.json').write_text(json.dumps(report,indent=2)+'\n')
    if rows:
        with (root/'trials.csv').open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    fetched=subprocess.run(['ssh','-S','/Users/xialeihuang/.ssh/torch-codex-cm',
        '-o','BatchMode=yes','xh2906@login.torch.hpc.nyu.edu',
        f'tar --ignore-failed-read -czf - -C {REMOTE} results selection.json submission.txt'],
        capture_output=True,check=True,timeout=180)
    with tarfile.open(fileobj=io.BytesIO(fetched.stdout)) as archive:
        for item in archive:
            if not item.isfile(): continue
            p=Path(item.name)
            assert not p.is_absolute() and '..' not in p.parts
            target=ROOT/p; target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(archive.extractfile(item).read())
    analyze(ROOT)
