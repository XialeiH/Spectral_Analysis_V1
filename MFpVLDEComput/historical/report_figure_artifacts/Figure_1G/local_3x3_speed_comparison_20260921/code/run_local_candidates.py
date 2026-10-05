"""Cold-process local DNN candidates; run only after the SNN comparison exits."""
import csv
import json
import os
import statistics
import subprocess
from pathlib import Path

ROOT = Path(os.environ.get('LOCAL_BENCH_ROOT', Path(__file__).resolve().parents[1]))
SOURCE = Path(os.environ.get('LOCAL_BENCH_SOURCE', ROOT.parent / 'tol5e3_3x3_4x4_20260921/local_smoke/source'))
OUT = ROOT / 'implementation_candidates'
MATLAB = '/Applications/MATLAB_R2025b.app/bin/matlab'
VARIANTS = ['reference', 'compact', 'paged', 'selected_paged_transposed']
REPEATS = [1, 4, 5]


def main():
    # Caller runs candidates serially after the preceding MATLAB process exits.
    assert (SOURCE / 'bundles/field_03HC.mat').is_file()
    env = os.environ.copy()
    env['SLURM_CPUS_PER_TASK'] = '8'  # Local thread-limit metadata, not a cluster allocation.
    for index, repeat in enumerate(REPEATS):
        variants = VARIANTS[index:] + VARIANTS[:index]
        for variant in variants:
            result = OUT / 'probe' / f'{variant}_repeat_{repeat:03d}.tsv'
            if result.exists():
                continue
            expression = (f"addpath('{OUT}/code'); run_impl_probe('{SOURCE}',"
                          f"'{OUT}',{repeat},'{variant}','probe');")
            print(f'START cold local {variant} repeat={repeat}', flush=True)
            with (OUT / 'logs' / f'{variant}_{repeat:03d}.log').open('w') as stream:
                subprocess.run([MATLAB, '-batch', expression], env=env,
                               stdout=stream, stderr=subprocess.STDOUT, check=True)
            with result.open() as stream:
                row = next(csv.DictReader(stream, delimiter='\t'))
            print(f"DONE {variant}:{row['Seconds']}s, passed={row['Passed']}", flush=True)
    summaries = {}
    for variant in VARIANTS:
        rows = []
        for repeat in REPEATS:
            with (OUT / 'probe' / f'{variant}_repeat_{repeat:03d}.tsv').open() as stream:
                rows.append(next(csv.DictReader(stream, delimiter='\t')))
        assert all(r['Passed'].lower() in ('1', 'true') and
                   int(r['Iterations']) == int(r['ReferenceIterations']) and
                   float(r['StepMaxAbsolute']) <= .001 and
                   float(r['StepMeanAbsolute']) <= .0001 and
                   float(r['StepResidualDifference']) <= .000001 for r in rows)
        times = [float(r['Seconds']) for r in rows]
        summaries[variant] = dict(n=len(times), mean=statistics.mean(times),
                                  sd=statistics.stdev(times), min=min(times), max=max(times),
                                  all_passed=True, individual_seconds=times,
                                  max_step_error=max(float(r['StepMaxAbsolute']) for r in rows))
    report = dict(platform='local Mac, eight threads', warmup=False, summaries=summaries,
                  under_point_one=[v for v in VARIANTS if summaries[v]['mean'] < .1],
                  note='Small cold-process screening study, not final100-repeat Torch data.')
    (OUT / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
