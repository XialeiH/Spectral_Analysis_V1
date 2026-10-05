"""Independent local smoke; never writes Torch result or figure directories."""
import csv
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path(os.environ.get('LOCAL_BENCH_ROOT', Path(__file__).resolve().parents[1]))
SOURCE = Path(os.environ.get('LOCAL_BENCH_SOURCE', ROOT.parent / 'tol5e3_3x3_4x4_20260921/local_smoke/source'))
MATLAB = '/Applications/MATLAB_R2025b.app/bin/matlab'
SHA = '19b51fb292c065ae59f0595aab9142bce49ab8ed342014de625a6eb490279516'
THREADS = 8


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matlab(expression, log):
    env = os.environ.copy()
    # The frozen rate-model wrapper reads this variable only as a thread limit.
    # This is a local process, not a Slurm allocation.
    env['SLURM_CPUS_PER_TASK'] = str(THREADS)
    with log.open('w') as stream:
        start = time.perf_counter()
        completed = subprocess.run([MATLAB, '-batch', expression], env=env,
                                   stdout=stream, stderr=subprocess.STDOUT, check=False)
        seconds = time.perf_counter() - start
    if completed.returncode:
        log.with_suffix('.failure.json').write_text(json.dumps(dict(
            exit_code=completed.returncode, failed_process_seconds=seconds,
            valid_benchmark=False), indent=2) + '\n')
        raise RuntimeError(f'MATLAB exit={completed.returncode}; see {log}')
    return seconds


def summarize():
    rows = []
    for path in sorted((ROOT / 'results').glob('*.tsv')):
        with path.open() as stream:
            rows.extend(csv.DictReader(stream, delimiter='\t'))
    summaries = {}
    for method in ['CG', 'DNN surrogate', 'SNN']:
        group = [r for r in rows if r['Method'] == method]
        if not group:
            continue
        values = [float(r['Seconds']) for r in group]
        sd = statistics.stdev(values) if len(values) > 1 else None
        summaries[method] = dict(n=len(values), mean=statistics.mean(values), sd=sd,
                                sem=sd / len(values)**.5 if sd is not None else None,
                                min=min(values), max=max(values), individual_seconds=values)
    result = dict(platform='local Mac, not Torch', threads=THREADS, summaries=summaries,
                  timing_scope={'CG': 'iteration and residual checks only',
                                'DNN surrogate': 'iteration and residual checks only',
                                'SNN': 'whole fresh MATLAB process, 15s simulation including I/O'})
    if 'CG' in summaries and 'DNN surrogate' in summaries:
        result['CG_over_DNN'] = summaries['CG']['mean'] / summaries['DNN surrogate']['mean']
    if 'SNN' in summaries and 'CG' in summaries:
        result['SNN_over_CG'] = summaries['SNN']['mean'] / summaries['CG']['mean']
    if 'SNN' in summaries and 'DNN surrogate' in summaries:
        result['SNN_over_DNN'] = summaries['SNN']['mean'] / summaries['DNN surrogate']['mean']
    result['complete'] = all(summaries.get(m, {}).get('n') == n
                             for m, n in [('CG', 5), ('DNN surrogate', 5), ('SNN', 1)])
    (ROOT / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


def main():
    driver = ROOT / 'code/Paper2_Fig7Comp_NW_LDE.m'
    assert sha(driver) == SHA
    inputs = ROOT / 'inputs'
    for name in ['AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat',
                 'LargeNWFixIni.mat', 'Paper2DriveNW_Conn.mat']:
        assert (inputs / name).is_file(), name
    assert (inputs / 'Utils/PoissonInputForNetwork.m').is_file()
    metadata = dict(platform=platform.platform(), machine=platform.machine(),
                    logical_cpus=os.cpu_count(), matlab=MATLAB, threads=THREADS,
                    created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    source=str(SOURCE), driver_sha256=SHA,
                    input_sha256={p.name: sha(p) for p in inputs.glob('*.mat')},
                    memory_bytes=int(subprocess.check_output(['sysctl', '-n', 'hw.memsize'])),
                    free_disk_bytes_before=shutil.disk_usage(ROOT).free,
                    note='No Slurm allocation; environment variable is reused only for MATLAB thread limit.')
    (ROOT / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for repeat in range(1, 6):
        methods = ['CG', 'DNN surrogate'] if repeat % 2 else ['DNN surrogate', 'CG']
        for method in methods:
            prefix = 'cg_dense' if method == 'CG' else 'dnn_fast'
            path = ROOT / 'results' / f'{prefix}_03HC_repeat_{repeat:03d}.tsv'
            if path.exists():
                continue
            print(f'START local {method} repeat={repeat}', flush=True)
            expression = (f"addpath('{ROOT}/code'); "
                          f"run_figure1g2_tol5e3('{SOURCE}','{ROOT}',3,{repeat},'{method}');")
            matlab(expression, ROOT / 'logs' / f'{prefix}_{repeat:03d}.log')
            with path.open() as stream:
                row = next(csv.DictReader(stream, delimiter='\t'))
            assert row['Converged'].lower() in ('1', 'true')
            assert float(row['FinalResidual']) < .005 and int(row['Iterations']) <= 200
            print(f"DONE local {method} repeat={repeat}: {row['Seconds']}s, {row['Iterations']} iterations", flush=True)
    summarize()

    result_path = ROOT / 'results/snn_15s_03HC_repeat_001.tsv'
    if result_path.exists():
        summarize()
        return
    # One SNN run writes ~5.5 GB; preserve a disk safety reserve.
    assert shutil.disk_usage(ROOT).free >= 8 * 1024**3, 'Insufficient free disk for full SNN trace; no SNN started.'
    trial = ROOT / 'snn_trial_001'
    data = trial / 'Data/Paper2_NetworkTuning/Fig1V4'
    save = data / 'Paper2NWSimulationData'
    for folder in [save, data / 'Paper2PlotingData', data / 'GlobConv', trial / 'Paper2Figs']:
        folder.mkdir(parents=True, exist_ok=True)
    if list(save.glob('DriveWkSp*.mat')):
        raise RuntimeError('Previous partial SNN output exists; inspect it rather than overwrite.')
    links = {trial / 'Utils': inputs / 'Utils',
             trial / 'LargeNWFixIni.mat': inputs / 'LargeNWFixIni.mat',
             data / 'AllMFPixPara_Paper2TuneFig1V4D2.mat': inputs / 'AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat',
             save / 'Paper2DriveNW_Conn.mat': inputs / 'Paper2DriveNW_Conn.mat'}
    for target, source in links.items():
        if not target.is_symlink():
            target.symlink_to(source)
    staged = trial / driver.name
    shutil.copy2(driver, staged)
    assert sha(staged) == SHA
    expression = (f"maxNumCompThreads({THREADS}); cd('{trial}'); run('{staged}'); "
                  "assert(N_HC==3 && T==15000 && dt==0.1 && TPar==15 && WinSize==9000); "
                  "assert(~BlowUp && TSec==15 && TimeN==10000); assert(maxNumCompThreads==8);")
    print('START local SNN:15s,150000updates,15 full segment saves', flush=True)
    seconds = matlab(expression, ROOT / 'logs/snn_15s_001.log')
    names = {f'DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_{i}s_NewSmear.mat' for i in range(1, 16)}
    actual = list(save.glob('DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat'))
    assert {p.name for p in actual} == names and all(p.stat().st_size > 0 for p in actual)
    assert (save / 'NWSimulationPix_0.0deg_NewSmear.mat').stat().st_size > 0
    assert sha(staged) == SHA
    row = dict(FieldHC=3, Repeat=1, Method='SNN', Seconds=seconds, SimulationMs=15000,
               dtMs=.1, Updates=150000, Threads=THREADS, Protocol='LocalOriginalPaper2Duration15s')
    with result_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    (ROOT / 'snn_validation.json').write_text(json.dumps(dict(
        valid=True, segments=len(actual), driver_sha256=sha(staged),
        total_segment_bytes=sum(p.stat().st_size for p in actual)), indent=2) + '\n')
    print(f'DONE local SNN:{seconds:.6f}s', flush=True)
    summarize()


if __name__ == '__main__':
    main()
