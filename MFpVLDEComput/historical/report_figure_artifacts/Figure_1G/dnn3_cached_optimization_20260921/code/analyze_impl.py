"""Summarize actual trials only; selection never uses target-derived times."""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

VARIANTS = ['paged', 'paged_quad', 'paged_padding', 'paged_aggregate', 'paged_combined']


def read_trials(root, phase):
    rows = []
    for path in sorted((root / phase).glob('*.tsv')):
        with path.open() as stream:
            entries = list(csv.DictReader(stream, delimiter='\t'))
        assert len(entries) == 1, str(path)
        rows.extend(entries)
    keys = [(r['Variant'], int(r['Repeat'])) for r in rows]
    assert len(keys) == len(set(keys)), 'Duplicate trial IDs'
    return rows


def valid(row):
    truth = lambda key: row[key].lower() in ('1', 'true')
    return (truth('Passed') and truth('Converged') and
            int(row['Iterations']) == int(row['ReferenceIterations']) and
            int(row['AllocatedCPUs']) == 16 and int(row['NNThreads']) == 8 and
            int(row['FieldHC']) == 3 and float(row['Seconds']) > 0 and
            float(row['MaxAbsolute']) <= 1e-3 and float(row['MeanAbsolute']) <= 1e-4 and
            float(row['StepMaxAbsolute']) <= 1e-3 and float(row['StepMeanAbsolute']) <= 1e-4 and
            float(row['ResidualDifference']) <= 1e-6 and float(row['StepResidualDifference']) <= 1e-6)


def summary(rows):
    times = [float(r['Seconds']) for r in rows]
    if not times:
        return {'n': 0}
    sd = statistics.stdev(times) if len(times) > 1 else None
    return dict(n=len(times), mean=statistics.mean(times), sd=sd,
                sem=sd / math.sqrt(len(times)) if sd is not None else None,
                min=min(times), max=max(times), all_passed=all(valid(r) for r in rows),
                max_step_error=max(float(r['StepMaxAbsolute']) for r in rows))


def complete(rows, count):
    return len(rows) == count and sorted(int(r['Repeat']) for r in rows) == list(range(1, count + 1))


def analyze(root, select=False):
    probe = read_trials(root, 'probe')
    groups = {v: [r for r in probe if r['Variant'] == v] for v in VARIANTS}
    probe_summary = {v: summary(rows) for v, rows in groups.items()}
    if select:
        assert all(complete(rows, 10) for rows in groups.values()), 'Probe incomplete'
        assert probe_summary['paged']['all_passed'], 'Invalid reference'
        eligible = [v for v in VARIANTS[1:] if probe_summary[v]['all_passed']]
        best = min(eligible, key=lambda v: probe_summary[v]['mean']) if eligible else 'none'
        if best != 'none' and probe_summary[best]['mean'] >= probe_summary['paged']['mean']:
            best = 'none'
        (root / 'selected_variant.txt').write_text(best + '\n')
    selected_path = root / 'selected_variant.txt'
    selected = selected_path.read_text().strip() if selected_path.exists() else None
    candidates = read_trials(root, 'validation')
    controls = read_trials(root, 'control')
    assert all(r['Variant'] == selected for r in candidates)
    candidate_summary = summary(candidates)
    control_summary = summary(controls)
    finished = complete(candidates, 100) and complete(controls, 100)
    ready = False
    if finished:
        ready = (candidate_summary['all_passed'] and control_summary['all_passed'] and
                 candidate_summary['mean'] < min(.37911415, .54671737, control_summary['mean']))
    result = dict(probe=probe_summary, selected_variant=selected, validation=candidate_summary,
                  control=control_summary, complete=finished, ready_to_plot=ready,
                  frozen_CG_mean=27.36131558, frozen_DNN4_mean=.54671737)
    if candidates:
        result['CG_over_DNN'] = 27.36131558 / candidate_summary['mean']
        result['target_100x_met'] = result['CG_over_DNN'] >= 100
    if candidates and controls:
        result['control_over_candidate_ratio_of_means'] = control_summary['mean'] / candidate_summary['mean']
    (root / 'implementation_status.json').write_text(json.dumps(result, indent=2) + '\n')
    for phase, rows in [('probe', probe), ('optimized', candidates), ('control', controls)]:
        if rows:
            with (root / (phase + '_trials.csv')).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('--select', action='store_true')
    args = parser.parse_args()
    analyze(args.root, args.select)
