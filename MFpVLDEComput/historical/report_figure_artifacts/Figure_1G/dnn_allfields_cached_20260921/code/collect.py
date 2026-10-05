"""Collect and validate full field groups; never infer or filter timing values."""
import csv
import io
import json
from pathlib import Path
import statistics
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
REMOTE='/scratch/xh2906/librarySCI_runs/figure1b_dnn_allfields_cached_20260921'
EXPECTED={4:100,6:20,8:20,10:20,20:20,30:20,40:20}


def valid(r,threads,variant):
    return (r['Passed'].lower() in ('1','true') and r['Converged'].lower() in ('1','true')
        and int(r['AllocatedCPUs'])==16 and int(r['NNThreads'])==threads
        and r['Variant']==variant and int(r['Iterations'])==int(r['ReferenceIterations'])
        and int(r['Iterations'])<=200 and float(r['FinalResidual'])<.005
        and float(r['Seconds'])>0 and float(r['MaxAbsolute'])<=.001
        and float(r['MeanAbsolute'])<=.0001 and float(r['StepMaxAbsolute'])<=.001
        and float(r['StepMeanAbsolute'])<=.0001 and float(r['ResidualDifference'])<=.000001
        and float(r['StepResidualDifference'])<=.000001)


def summary(rows):
    if not rows: return {'n':0}
    values=[float(r['Seconds']) for r in rows]
    sd=statistics.stdev(values) if len(values)>1 else None
    return dict(n=len(values),mean=statistics.mean(values),sd=sd,
        sem=sd/len(values)**.5 if sd is not None else None,min=min(values),max=max(values))


def read(root,phase):
    rows=[]
    for p in sorted((root/phase).glob('*.tsv')):
        with p.open() as f: data=list(csv.DictReader(f,delimiter='\t'))
        assert len(data)==1
        rows.extend(data)
    keys=[(r['FieldHC'],r['Repeat']) for r in rows]
    assert len(keys)==len(set(keys))
    return rows


def analyze(root):
    gate=read(root,'gate'); candidates=read(root,'validation'); controls=read(root,'control')
    fields={}; errors=[]
    for field,n in EXPECTED.items():
        g=[r for r in gate if int(r['FieldHC'])==field]
        a=sorted([r for r in candidates if int(r['FieldHC'])==field],key=lambda r:int(r['Repeat']))
        b=sorted([r for r in controls if int(r['FieldHC'])==field],key=lambda r:int(r['Repeat']))
        checks=all(valid(r,8,'paged_combined') for r in g+a) and all(valid(r,16,'reference') for r in b)
        paired=True
        for c in a:
            matches=[r for r in b if r['Repeat']==c['Repeat']]
            if matches and matches[0]['Node']!=c['Node']: paired=False
        finished=(sorted(int(r['Repeat']) for r in g)==[1,4]
            and [int(r['Repeat']) for r in a]==list(range(1,n+1))
            and [int(r['Repeat']) for r in b]==list(range(1,n+1)))
        if not checks or not paired: errors.append(f'{field}HC validation failed')
        report=dict(candidate=summary(a),control=summary(b),gate_n=len(g),
            checks_passed=checks and paired,complete=finished,ready=finished and checks and paired)
        if a and b: report['control_over_candidate']=report['control']['mean']/report['candidate']['mean']
        fields[str(field)]=report
    result=dict(fields=fields,errors=errors,complete=all(v['complete'] for v in fields.values()),
        ready_to_plot=all(v['ready'] for v in fields.values()),
        note='Publish all complete validated field groups without filtering slower timings.')
    (root/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,rows in [('optimized_trials',candidates),('control_trials',controls),('gate_trials',gate)]:
        if rows:
            with (root/(name+'.csv')).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    script="""
import pathlib,sys,tarfile
r=pathlib.Path(%r)
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|') as t:
    for phase in ['gate','validation','control']:
        for p in sorted((r/phase).glob('*.tsv')): t.add(p,arcname=str(p.relative_to(r)))
    if (r/'submission.txt').exists(): t.add(r/'submission.txt',arcname='submission.txt')
""" % REMOTE
    command="bash -lc 'module load anaconda3/2025.06; source /share/apps/anaconda3/2025.06/etc/profile.d/conda.sh; conda activate py310; python -'"
    fetched=subprocess.run(['ssh','-S','/Users/xialeihuang/.ssh/torch-codex-cm','-o','BatchMode=yes',
        'xh2906@login.torch.hpc.nyu.edu',command],input=script.encode(),capture_output=True,check=True,timeout=180)
    with tarfile.open(fileobj=io.BytesIO(fetched.stdout)) as archive:
        for item in archive:
            p=Path(item.name)
            assert item.isfile() and not p.is_absolute() and '..' not in p.parts
            target=ROOT/p; target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(archive.extractfile(item).read())
    analyze(ROOT)
