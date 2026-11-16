"""Prepare supporting pulse diagnostics from existing saved trajectories.

Run ``python3 Gain_Switching/scripts/sm_figures.py`` from any directory.
Requires NumPy and gnuplot. Metrics use every saved sample; only the plotting
exports are decimated. Original data and previously generated figures are kept.
The normalization and pulse extraction reuse gain_dependence.py.
"""
from pathlib import Path
import csv
import json
import subprocess

import numpy as np
import gain_dependence as gd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Gain_Switching/Pub_Ready'
CASES = [
    (1.15, ROOT / 'NEW_stuff/Data/GainDependence/1.15Gth.dat'),
    (1.20, ROOT / 'NEW_stuff/Extra/1.2Gth.dat'),
    (1.25, ROOT / 'NEW_stuff/Extra/1.25Gth.dat'),
]
PERIOD = gd.PERIOD_PS
DISPLAY_STRIDE = 30


def crossings(x, y, level, rising):
    """Linearly interpolated crossings, preserving the exported time precision."""
    hit = np.flatnonzero((y[:-1] < level) & (y[1:] >= level) if rising
                         else (y[:-1] >= level) & (y[1:] < level))
    return [float(x[i] + (level-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))
            for i in hit]


def main():
    OUT.mkdir(exist_ok=True)
    (ROOT / 'LaTeX').mkdir(exist_ok=True)
    rows, records, summary = [], [], []
    gth = float(np.median([gd.header_gain(path)/gain for gain, path in CASES]))
    for gain, path in CASES:
        if abs(gd.header_gain(path)/gth-gain) > .0001:
            raise ValueError(f'Gain header inconsistent with filename: {path}')
        t, modulus, n = np.loadtxt(path, usecols=(0, 5, 6), unpack=True)
        if not np.isfinite([t, modulus, n]).all() or np.any(np.diff(t) < 0) or t[-1] < 200:
            raise ValueError(f'Invalid/incomplete trajectory: {path}')
        intensity = (modulus/gd.esat(gain))**2
        case_rows = []
        for cycle in range(4):
            mask = (t >= PERIOD*cycle) & (t < PERIOD*(cycle+1))
            x, y, nn = t[mask]-PERIOD*cycle, intensity[mask], n[mask]
            metric = gd.analyze_cycle(t, intensity, PERIOD*cycle)
            if not metric['clear_pulse']:
                raise ValueError(f'No clear pulse at gain {gain}, cycle {cycle+1}')
            peak = int(np.argmax(y))
            rise = [v for v in crossings(x, y, y.max()/2, True) if v <= x[peak]][-1]
            fall = [v for v in crossings(x, y, y.max()/2, False) if v >= x[peak]][0]
            if not np.isclose(fall-rise, metric['fwhm_ps'], atol=1e-12):
                raise ValueError('Half-maximum extraction differs from main figure')
            n_up, n_down = crossings(x, nn, 1/gain, True), crossings(x, nn, 1/gain, False)
            row = dict(gain_ratio=gain, cycle=cycle+1,
                       peak_time_ps=metric['peak_time_ps'], FWHM_ps=metric['fwhm_ps'],
                       peak_Iem_over_Isat=metric['peak_intensity'],
                       half_rise_ps=rise, half_fall_ps=fall,
                       N_at_peak=float(nn[peak]), Ncrit=1/gain,
                       N_up_ps=n_up[0] if n_up else None,
                       N_down_ps=n_down[0] if n_down else None,
                       N_crossing_count=len(n_up)+len(n_down))
            rows.append(row)
            case_rows.append(row)
            # Include the exact stored peak and nearest crossing samples in the
            # visual export, in addition to regularly spaced display samples.
            landmarks = [rise, fall, 25] + n_up + n_down
            near = [int(np.argmin(abs(x-v))) for v in landmarks]
            indices = np.unique(np.r_[np.arange(0, len(x), DISPLAY_STRIDE),
                                       peak, len(x)-1, near])
            target = .5*(1-np.cos(2*np.pi*x/PERIOD))
            np.savetxt(OUT/f'sm_transition_{gain:.2f}_cycle{cycle+1}.dat',
                       np.column_stack((x[indices], nn[indices],
                                        target[indices], 1e3*y[indices])),
                       fmt='%.12g', header='time_within_cycle_ps N Ntarget 1000*Iem/Isat')
        records.append(dict(gain_ratio=gain, source=str(path.relative_to(ROOT)),
                            sha256=gd.digest(path), Esat=float(gd.esat(gain)),
                            total_saved_samples=len(t)))
        last = case_rows[-1]
        summary.append(dict(gain_ratio=gain,
                            final_peak_time_ps=last['peak_time_ps'],
                            final_FWHM_ps=last['FWHM_ps'],
                            final_peak_Iem_over_Isat=last['peak_Iem_over_Isat'],
                            final_half_rise_ps=last['half_rise_ps'],
                            final_half_fall_ps=last['half_fall_ps'],
                            cycle_peak_time_span_ps=float(np.ptp([r['peak_time_ps'] for r in case_rows])),
                            cycle_FWHM_span_ps=float(np.ptp([r['FWHM_ps'] for r in case_rows])),
                            cycle_peak_intensity_relative_span_percent=100*float(np.ptp([r['peak_Iem_over_Isat'] for r in case_rows]))/last['peak_Iem_over_Isat']))
    for name, table in [('sm_transition_metrics.csv', rows), ('sm_transition_summary.csv', summary)]:
        with (OUT/name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
    # Gnuplot includes the final-cycle marker coordinates, from full-sample metrics.
    with (OUT/'sm_transition_parameters.gp').open('w') as stream:
        for name, field, multiplier in [('peak_times', 'peak_time_ps', 1),
                                        ('peaks', 'peak_Iem_over_Isat', 1000),
                                        ('half_rises', 'half_rise_ps', 1),
                                        ('half_falls', 'half_fall_ps', 1)]:
            stream.write(f"array {name}[3] = [" + ', '.join(f'{r[field]*multiplier:.12g}' for r in rows if r['cycle'] == 4) + ']\n')
    (OUT/'sm_transition_provenance.json').write_text(json.dumps(dict(
        files=records, period_ps=PERIOD, display_window_ps=[12, 38],
        mechanism_cycle_ps=[150, 200], overlay_cycles=[[0, 50], [50, 100], [100, 150], [150, 200]],
        normalization='Esat=0.317931*sqrt(1.10/g); same gain_dependence.esat as the main figure',
        extraction='gain_dependence.analyze_cycle on all stored samples; linearly interpolated half-maximum crossings bracketing each maximum',
        display=f'Every {DISPLAY_STRIDE}th stored sample plus maxima, endpoints, and nearest crossing samples; plotting decimation only',
        pump='Ntarget=0.5*(1-cos(2*pi*t/50))',
        note='Cycle agreement measures repeatability in exported deterministic trajectories, not time-step convergence, stochastic timing uncertainty, or stability under arbitrary perturbations. Rounded timestamps limit the timing precision. No simulation is run.',
        plot_script='Gain_Switching/gnuplot/sm_transition.gp',
        outputs=['LaTeX/sm_gain_switching_transition_mechanism.pdf', 'LaTeX/sm_gain_switching_transition_cycles.pdf'],
        summary=summary), indent=2) + '\n')
    subprocess.run(['gnuplot', str(ROOT/'Gain_Switching/gnuplot/sm_transition.gp')], check=True)
    for row in summary:
        print(row)


if __name__ == '__main__':
    main()
