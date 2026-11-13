#!/usr/bin/env python3
"""Summarize CUP's recorded decisions; no independent convergence detector."""
import csv
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paper = root / 'data/output/[nanoph-2024-0491]_rangeE0'
work = root / 'data/output/convergence_validation'
isat = float(next(line.split('=')[1] for line in (work/'inputs/parameters.gp').read_text().splitlines() if line.startswith('Isat =')))
checks = []
for seed in ('08', '09', '10', '11'):
    original = (paper/f'numtime{seed}.dat').read_bytes()
    rebuilt = (work/f'fixed/{seed}/numtime.dat').read_bytes()
    checks.append(f'{seed}\t{original == rebuilt}\t{hashlib.sha256(original).hexdigest()}\t{hashlib.sha256(rebuilt).hexdigest()}')
(work/'fixed_comparison.tsv').write_text('seed\tbyte_identical\tbaseline_sha256\trebuilt_sha256\n'+'\n'.join(checks)+'\n')
if not all(line.split('\t')[1] == 'True' for line in checks):
    raise SystemExit('Fixed-duration mismatch: inspect fixed_comparison.tsv')
rows = []
for case in ('automatic', 'tight'):
    for seed in ('08', '09', '10', '11'):
        folder = work/case/seed
        metadata, samples = {}, []
        for line in (folder/'numtime.dat').read_text().splitlines():
            if line.startswith('#'):
                if ' = ' in line:
                    key, value = line[2:].split(' = ', 1)
                    metadata[key] = value
            elif line.strip():
                values = line.split()
                t, re, im = map(float, (values[0], values[3], values[4]))
                samples.append((t, re*re+im*im))
        state = metadata['termination_reason']
        if state not in ('sustained_flat', 'sustained_periodic'):
            raise SystemExit(f'{case}/{seed} did not converge: {state}')
        start, end, mean, pp = (float(metadata[f'convergence_{key}']) for key in ('start_ps', 'confirmation_end_ps', 'mean_intensity', 'peak_to_peak'))
        rows.append([case, f'1e-{int(seed)}', state, start, float(metadata['convergence_qualifying_end_ps']), end,
                     float(metadata['convergence_detected_at_ps']), mean/isat, pp/isat, 100*pp/mean,
                     float(metadata['convergence_mean_period_ps'])])
        # Saved output has 0.1 ps cadence and rounded field columns. The summary
        # above uses full-integration-resolution detector statistics, not these rows.
        with (folder/'accepted_tail.dat').open('w') as out:
            out.write('# time_ps relative_intensity_deviation_percent; accepted+confirmation interval\n')
            for t, intensity in samples:
                if start <= t <= end:
                    out.write(f'{t:.17g} {(intensity/mean-1)*100:.12g}\n')
with (work/'summary.tsv').open('w') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['case', 'raw_E0', 'state', 'start_ps', 'qualifying_end_ps', 'confirmation_end_ps', 'stop_ps',
                     'mean_I_over_Isat', 'peak_to_peak_I_over_Isat', 'relative_peak_to_peak_percent', 'period_ps'])
    writer.writerows(rows)
print('All four fixed-duration traces are byte-identical to the saved reproduction.')
for row in rows:
    print(f'{row[0]:9s} E0={row[1]:5s}: stop={row[6]:.3f} ps; p-p/mean={row[9]:.6f}%; period={row[10]:.6f} ps')
