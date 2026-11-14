"""Recompute pulse metrics and render gain_dependence.{pdf,png}.

Run with python3 from any working directory. Requires numpy and gnuplot.
Input files stay in NEW_stuff/Data/GainDependence; no simulation is rerun.
The filename supplies the intended gain; headers are checked for consistency.
Duplicate files and conflicting gain headers are flagged, not silently omitted.
Normalization uses the simulation-unit formula Esat=1/sqrt(abs(G)*tau1)
documented in Data_Description.pdf, pp. 6-7, at fixed tau1. The rounded
Esat=0.317931 at G/Gth=1.10 anchors its numerical scale. This convention
is not an assertion that the physical saturation field varies with gain.
"""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess
import numpy as np

BASE = Path(__file__).resolve().parents[1]
INPUT = BASE.parent / 'NEW_stuff/Data/GainDependence'
OUTPUT = BASE / 'Pub_Ready'
PERIOD_PS = 50.0
STABILITY_PERCENT = 1.0
ESAT_DOCUMENTED = {1.05: 0.325413, 1.10: 0.317931, 1.15: 0.310942}


def esat(gain):
    """Documented fixed-tau1 law in simulation units, anchored at g=1.10."""
    if not np.isfinite(gain) or gain <= 0:
        raise ValueError('Gain ratio must be finite and positive')
    return 0.317931 * np.sqrt(1.10 / gain)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def header_gain(path):
    with path.open() as stream:
        for line in stream:
            if not line.startswith('#'):
                break
            match = re.match(r'#\s*G\s*=\s*([\d.eE+-]+)', line)
            if match:
                return float(match[1])
    return float('nan')


def analyze_cycle(t, intensity, start):
    mask = (t >= start) & (t < start + PERIOD_PS)
    x, y = t[mask], intensity[mask]
    if len(x) < 3:
        raise ValueError('Insufficient samples in a complete cycle')
    peak = int(np.argmax(y))
    result = dict(cycle_start_ps=float(start), peak_intensity=float(y[peak]),
                  peak_time_ps=float(x[peak] - start), fwhm_ps=float('nan'),
                  clear_pulse=False, reason='')
    if peak in (0, len(y)-1) or y[peak] <= 0:
        result['reason'] = 'Nonpositive peak or peak at cycle boundary'
        return result
    if (y[peak] - np.percentile(y, 10)) / y[peak] < 0.10:
        result['reason'] = 'Peak-background contrast below 10 percent'
        return result
    half = y[peak] / 2
    idx = np.arange(len(y)-1)
    left = np.flatnonzero((y[:-1] < half) & (y[1:] >= half) & (idx < peak))
    right = np.flatnonzero((y[:-1] >= half) & (y[1:] < half) & (idx >= peak))
    if not len(left) or not len(right):
        result['reason'] = 'Missing half-maximum crossing'
        return result
    def crossing(i):
        return x[i] + (half-y[i]) * (x[i+1]-x[i]) / (y[i+1]-y[i])
    width = float(crossing(right[0]) - crossing(left[-1]))
    result.update(fwhm_ps=width, clear_pulse=width > 0,
                  reason='' if width > 0 else 'Nonpositive width')
    return result


def difference(a, b):
    denominator = (abs(a) + abs(b))/2
    return float(100*abs(a-b)/denominator) if denominator else 0.0


def main():
    # All three supplied anchors must agree within their six-decimal rounding.
    for gain, documented in ESAT_DOCUMENTED.items():
        if not np.isclose(esat(gain), documented, rtol=3e-6, atol=0):
            raise ValueError(f'Normalization conflicts with documented Esat at {gain}')
    OUTPUT.mkdir(exist_ok=True)
    files = []
    for path in INPUT.glob('*.dat'):
        match = re.fullmatch(r'([\d.]+)Gth', path.stem)
        if match:
            files.append((float(match[1]), path))
    files.sort()
    expected = [round(1.05+i*0.01, 2) for i in range(11)]
    if [g for g, _ in files] != expected:
        raise ValueError('Expected exactly the eleven gain files from 1.05 to 1.15')
    # Infer the common threshold robustly from rounded gain headers.
    gth = float(np.nanmedian([header_gain(p)/g for g, p in files]))
    seen, seen_q, rows, provenance = {}, {}, [], []
    for gain, path in files:
        warnings = []
        sha = digest(path)
        if sha in seen:
            warnings.append('Duplicate of ' + seen[sha])
        else:
            seen[sha] = path.name
        qpath = path.with_name(path.stem + 'q.dat')
        qsha = digest(qpath) if qpath.exists() else None
        if qsha is not None:
            if qsha in seen_q:
                warnings.append('Companion duplicate of ' + seen_q[qsha])
            else:
                seen_q[qsha] = qpath.name
        actual_gain = header_gain(path)
        if not np.isfinite(actual_gain) or abs(actual_gain/gth-gain) > 0.0001:
            warnings.append('Gain header conflicts with filename')
        t, modulus = np.loadtxt(path, usecols=(0, 5), unpack=True)
        if not np.isfinite(t).all() or not np.isfinite(modulus).all() or np.any(np.diff(t) < 0):
            raise ValueError(f'Invalid samples in {path}')
        # Anchor cycles to pump t=0, not the first saved sample. Include the
        # last complete cycle; exclude its shared right endpoint in metrics.
        end = np.floor(t[-1]/PERIOD_PS)*PERIOD_PS
        if end < 2*PERIOD_PS:
            raise ValueError(f'Fewer than two complete cycles in {path}')
        saturation = float(esat(gain))
        intensity = (modulus/saturation)**2
        previous = analyze_cycle(t, intensity, end-2*PERIOD_PS)
        last = analyze_cycle(t, intensity, end-PERIOD_PS)
        changes = {key: difference(previous[key], last[key])
                   for key in ('peak_intensity', 'fwhm_ps', 'peak_time_ps')}
        stable = bool(previous['clear_pulse'] and last['clear_pulse'] and
                      all(value <= STABILITY_PERCENT for value in changes.values()))
        if not stable:
            warnings.append('No stable clear pulse in last two cycles')
        row = dict(gain_ratio=gain, Esat=saturation,
                   normalization='documented fixed-tau1 scaling',
                   peak_Iem_over_Isat=last['peak_intensity'], fwhm_ps=last['fwhm_ps'],
                   peak_time_in_cycle_ps=last['peak_time_ps'],
                   cycle_start_ps=last['cycle_start_ps'], cycle_end_ps=float(end),
                   clear_pulse=last['clear_pulse'], pulse_reason=last['reason'],
                   peak_change_percent=changes['peak_intensity'],
                   fwhm_change_percent=changes['fwhm_ps'],
                   peak_time_change_percent=changes['peak_time_ps'], stable=stable,
                   unverified=bool(warnings), warning='; '.join(warnings),
                   source=path.name, header_G=actual_gain, sha256=sha)
        rows.append(row)
        provenance.append(dict(source=str(path.relative_to(BASE.parent)),
                               companion_sha256=qsha, previous_cycle=previous, last_cycle=last))
        print(f'{gain:.2f}: FWHM={last["fwhm_ps"]:.6f} ps, '
              f'peak={last["peak_intensity"]:.8g}, time={last["peak_time_ps"]:.3f} ps'
              + (' — WARNING: ' + row['warning'] if warnings else ''))
    with (OUTPUT/'gain_dependence_summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (OUTPUT/'gain_dependence_plot.dat').open('w') as stream:
        stream.write('# gain FWHM_ps peak_Iem_over_Isat peak_time_ps unverified\n')
        for row in rows:
            stream.write(f'{row["gain_ratio"]:.2f} {row["fwhm_ps"]:.12g} '
                         f'{row["peak_Iem_over_Isat"]:.12g} '
                         f'{row["peak_time_in_cycle_ps"]:.12g} {int(row["unverified"])}\n')
    (OUTPUT/'gain_dependence_provenance.json').write_text(json.dumps(dict(
        period_ps=PERIOD_PS, stability_tolerance_percent=STABILITY_PERCENT,
        normalization='All gains: Esat=0.317931*sqrt(1.10/g), g=Gmax/Gth',
        normalization_source='Data_Description.pdf, pp. 6-7: Esat=1/sqrt(abs(G)*tau1)',
        normalization_convention='Simulation units; fixed tau1; rounded anchor at g=1.10',
        normalization_anchors=ESAT_DOCUMENTED,
        physical_units_note=('doc/ns_dynamic.pdf Eq. (12) and p. 45 define the physical '
                             'saturation field. Conversion from simulation units remains '
                             'to be checked against the simulation source.'),
        inferred_threshold_gain=gth, files=provenance), indent=2))
    subprocess.run(['gnuplot', str(BASE/'gnuplot/gain_dependence.gp')], check=True)


if __name__ == '__main__':
    main()
