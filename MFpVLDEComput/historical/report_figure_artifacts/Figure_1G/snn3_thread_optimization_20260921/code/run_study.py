"""Unmodified SNN20 driver, fresh-process timing, isolated thread-only study."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829')
INPUTS = BASE / 'snn_paper2_original_exact_20260901_160003/inputs'
SHA = 'b3958be55c41dba44f674de73e3dfb8bf37e2b552c7f4b145b824872197f4c6b'
THREADS = [1, 2, 4, 8, 16]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trial(phase, threads, repeat):
    name = f'{phase}_t{threads:02d}_r{repeat:03d}'
    work = ROOT / 'trials' / name
    work.mkdir(parents=True, exist_ok=False)
    data = work / 'Data/Paper2_NetworkTuning/Fig1V4'
    save = data / 'Paper2NWSimulationData'
    for p in [save, data/'Paper2PlotingData', data/'GlobConv', work/'Paper2Figs']:
        p.mkdir(parents=True, exist_ok=True)
    links = {work/'Utils': BASE/'snn_repo/Utils',
             work/'LargeNWFixIni.mat': INPUTS/'LargeNWFixIni.mat',
             data/'AllMFPixPara_Paper2TuneFig1V4D2.mat': INPUTS/'AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat',
             save/'Paper2DriveNW_Conn.mat': INPUTS/'Paper2DriveNW_Conn.mat'}
    for target, source in links.items():
        assert source.exists()
        target.symlink_to(source)
    driver = work/'Paper2_Fig7Comp_NW_LDE.m'
    shutil.copy2(ROOT/'code'/driver.name, driver)
    assert digest(driver) == SHA
    # Preserve the default fresh-MATLAB random sequence for every configuration.
    expr = (f"maxNumCompThreads({threads}); rng('default'); cd('{work}'); "
            f"run('{driver}'); assert(N_HC==3 && T==20000 && dt==.1 && TPar==20 && WinSize==9000); "
            f"assert(~BlowUp && TSec==20 && TimeN==10000 && maxNumCompThreads=={threads});")
    with (ROOT/'logs'/f'{name}.matlab.log').open('w') as log:
        start = time.perf_counter()
        process = subprocess.run(['matlab','-batch',expr], stdout=log, stderr=subprocess.STDOUT)
        seconds = time.perf_counter()-start
    if process.returncode:
        (ROOT/'results'/f'{name}.failure.json').write_text(json.dumps(dict(
            exit_code=process.returncode, failed_seconds=seconds, valid=False)))
        raise RuntimeError(f'MATLAB failed: {name}')
    expected = {f'DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_{i}s_NewSmear.mat' for i in range(1,21)}
    files = list(save.glob('DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat'))
    assert {p.name for p in files} == expected and all(p.stat().st_size>0 for p in files)
    assert (save/'NWSimulationPix_0.0deg_NewSmear.mat').stat().st_size>0
    assert digest(driver)==SHA and int(os.environ['SLURM_CPUS_PER_TASK'])==16
    row = dict(Phase=phase,Threads=threads,Repeat=repeat,Seconds=seconds,FieldHC=3,
               SimulationMs=20000,dtMs=.1,Updates=200000,Segments=20,CPUs=16,
               Seed='default',Node=os.environ.get('SLURMD_NODENAME',''),DriverSHA=SHA,Trial=name)
    with (ROOT/'results'/f'{name}.tsv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(row),delimiter='\t')
        writer.writeheader(); writer.writerow(row)
    print(json.dumps(row),flush=True)
    return work


def rows():
    result=[]
    for p in sorted((ROOT/'results').glob('*.tsv')):
        with p.open() as f:
            data=list(csv.DictReader(f,delimiter='\t'))
        assert len(data)==1
        result.extend(data)
    return result


def select():
    probe=[r for r in rows() if r['Phase']=='probe']
    means={}
    for threads in THREADS:
        group=[r for r in probe if int(r['Threads'])==threads]
        assert sorted(int(r['Repeat']) for r in group)==[1,2,3]
        means[threads]=statistics.mean(float(r['Seconds']) for r in group)
    best=min(means,key=means.get)
    (ROOT/'selection.json').write_text(json.dumps(dict(threads=best,means=means),indent=2))
    print('Selected',best,'threads',means,flush=True)


def validate(repeat):
    selected=json.loads((ROOT/'selection.json').read_text())['threads']
    if selected==16:
        print('No faster probe configuration than16 threads; no validation runs needed.')
        return
    configs=[('candidate',selected),('control',16)]
    if repeat%2==0: configs.reverse()
    completed={phase:trial(phase,threads,repeat) for phase,threads in configs}
    expr=(f"addpath('{ROOT}/code'); validate_outputs('{completed['candidate']}',"
          f"'{completed['control']}','{ROOT}/results/pair_{repeat:03d}.json');")
    with (ROOT/'logs'/f'pair_{repeat:03d}.log').open('w') as log:
        subprocess.run(['matlab','-batch',expr],stdout=log,stderr=subprocess.STDOUT,check=True)


if __name__=='__main__':
    phase=sys.argv[1]
    if phase=='select': select()
    elif phase=='validate': validate(int(os.environ['SLURM_ARRAY_TASK_ID']))
    elif phase=='probe':
        index=int(os.environ['SLURM_ARRAY_TASK_ID'])-1
        repeat=index//5+1
        threads=THREADS[(index%5+repeat-1)%5]
        trial('probe',threads,repeat)
    else: raise ValueError(phase)
