#!/usr/bin/env python3
"""Summarize measured CP timing trials and check paired M comparisons."""
from pathlib import Path
from collections import defaultdict
import csv, math
ROOT = Path(__file__).resolve().parents[1]

def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))

def write(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def summarize_sync(phase):
    groups = defaultdict(list)
    for r in read(ROOT / 'results' / f'sync_{phase}_trials.csv'):
        groups[(r['modulation'], int(r['snr_db']), int(r['symbols']))].append(r)
    rows = []
    for (mode, snr, m), trial in sorted(groups.items()):
        assert len(trial) == 2000 and {int(r['trial']) for r in trial} == set(range(2000))
        assert all(sum(int(r['true_cp']) == d for r in trial) == 25 for d in range(80))
        errors = [int(r['estimated_cp']) - int(r['true_cp']) for r in trial]
        assert all(e == int(r['error']) for e, r in zip(errors, trial))
        assert all(int(r['correct']) == (e == 0) for e, r in zip(errors, trial))
        circular = [int(r['circular_error']) for r in trial]
        good = errors.count(0)
        rows.append(dict(phase=phase, modulation=mode, snr_db=snr, symbols=m, trials=2000,
                         correct=good, failures=2000-good, correct_rate=good/2000,
                         failure_rate=1-good/2000, rmse=math.sqrt(sum(e*e for e in errors)/2000),
                         mae=sum(abs(e) for e in errors)/2000,
                         circular_rmse=math.sqrt(sum(e*e for e in circular)/2000),
                         required_samples=80*m+79, direct_pairs=80*16*m))
    assert len(rows) == (66 if phase == 'selection' else 198)
    write(ROOT / 'results' / f'sync_{phase}_summary.csv', rows)
    print(f'{phase}: {sum(r["trials"] for r in rows):,} estimates, {len(rows)} conditions')

def summarize_multipath():
    groups = defaultdict(list)
    for r in read(ROOT / 'results/m40/multipath_timing_trials.csv'):
        groups[(r['modulation'], r['profile'], int(r['snr_db']), int(r['symbols']))].append(r)
    rows = []
    for (mode, profile, snr, m), rs in sorted(groups.items()):
        assert len(rs) == 1024
        e = [int(r['delta']) for r in rs]
        rows.append(dict(modulation=mode, profile=profile, snr_db=snr, symbols=m, trials=1024,
                         exact_cp=e.count(0), valid_fft=sum(int(r['fft_window_valid']) for r in rs),
                         mean_shift=sum(e)/1024, rmse=math.sqrt(sum(x*x for x in e)/1024)))
    write(ROOT / 'results/m40/multipath_timing_summary.csv', rows)
    print('multipath: 24 paired conditions')

def check_ber():
    b64 = read(ROOT / 'results/ber.csv')
    b40 = read(ROOT / 'results/m40/ber.csv')
    assert len(b64) == len(b40) == 18
    for a, b in zip(b64, b40):
        assert a['modulation'] == b['modulation'] and a['snr_db'] == b['snr_db']
        assert a['errors'] == b['errors'] and a['bits'] == b['bits'] == '33554432'
        for r in (a, b):
            assert r['estimated_cp'] == r['true_cp']
            assert abs(float(r['ber'])-int(r['errors'])/int(r['bits'])) < 1e-15
    print('AWGN BER: 18 paired points with equal counts')

if __name__ == '__main__':
    for phase in ('selection', 'fine'):
        summarize_sync(phase)
    summarize_multipath()
    check_ber()
