"""Recompute the original, extended, and additional gain series and render the linear figure."""
from pathlib import Path
import csv
import json
import subprocess
import numpy as np
import gain_dependence as gd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parent
OUT = BASE / 'Pub_Ready'
STEM = 'gain_dependence_linear_extended'

def main():
    files = []
    for folder in (ROOT/'NEW_stuff/Data/GainDependence', ROOT/'NEW_stuff/Extra', ROOT/'EXTRA_stuff'):
        for path in folder.glob('*Gth.dat'):
            files.append((float(path.stem[:-3]), path))
    files.sort()
    expected = [round(1.05+i*.01, 2) for i in range(11)] + [round(1.20+i*.05, 2) for i in range(7)]
    expected = sorted(expected + [1.16, 1.17, 1.18, 1.19, 1.21, 1.22, 1.23, 1.24, 1.26, 1.27, 1.28, 1.29, 1.33, 1.38, 1.43, 1.48])
    if [g for g, p in files] != expected:
        raise ValueError('Unexpected gain coverage')
    gth = np.median([gd.header_gain(p)/g for g,p in files])
    rows, seen = [], set()
    for gain, path in files:
        sha = gd.digest(path)
        if sha in seen or not np.isfinite(gd.header_gain(path)) or abs(gd.header_gain(path)/gth-gain) > .0001:
            raise ValueError(f'Duplicate or mismatched gain: {path}')
        seen.add(sha)
        # The exported dipole coefficient equals alpha times the applied field.
        # Check the normalized drive from the first saved sample and the
        # resonant energy recorded in the header without assuming filenames
        # alone identify the run configuration.
        with path.open() as stream:
            header = []
            for line in stream:
                if line.startswith('#'):
                    header.append(line)
                elif line.strip():
                    first = np.fromstring(line, sep=' ')
                    break
        energies = [float(line.split('=')[1].split()[0]) for line in header
                    if line.startswith('# ome =') or line.startswith('# ome_em =')]
        drive = abs(complex(first[3], first[4]) / complex(first[1], first[2])) / gd.esat(gain)
        if len(energies) != 2 or not np.allclose(energies, 2.8122, rtol=0, atol=1e-5):
            raise ValueError(f'Nonresonant or missing frequency header: {path}')
        if not np.isclose(drive, 1e-9, rtol=1e-4, atol=0):
            raise ValueError(f'Unexpected normalized probe amplitude: {path}: {drive}')
        t, modulus = np.loadtxt(path, usecols=(0,5), unpack=True)
        if not np.isfinite(t).all() or not np.isfinite(modulus).all() or np.any(np.diff(t)<0) or t[0] > 100 or t[-1] < 200 or np.any(modulus < 0):
            raise ValueError(f'Invalid trajectory: {path}')
        intensity = (modulus/gd.esat(gain))**2
        previous = gd.analyze_cycle(t,intensity,100)
        last = gd.analyze_cycle(t,intensity,150)
        changes = {k:gd.difference(previous[k],last[k]) for k in ('fwhm_ps','peak_intensity','peak_time_ps')}
        if not previous['clear_pulse'] or not last['clear_pulse'] or max(changes.values()) > gd.STABILITY_PERCENT:
            raise ValueError(f'Unstable pulse: {path}')
        rows.append(dict(gain_ratio=gain, fwhm_ps=last['fwhm_ps'], peak_Iem_over_Isat=last['peak_intensity'], peak_time_ps=last['peak_time_ps'], Esat=float(gd.esat(gain)), source=str(path.relative_to(ROOT)), sha256=sha, changes_percent=changes))
    with (OUT/(STEM+'_plot.dat')).open('w') as f:
        f.write('# gain FWHM_ps peak_Iem_over_Isat peak_time_ps unverified\n')
        for r in rows:
            f.write(f"{r['gain_ratio']:.2f} {r['fwhm_ps']:.12g} {r['peak_Iem_over_Isat']:.12g} {r['peak_time_ps']:.12g} 0\n")
    with (OUT/(STEM+'_summary.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (OUT/(STEM+'_provenance.json')).write_text(json.dumps(dict(cycle_ps=[150,200], normalization='Esat=0.317931*sqrt(1.10/gain_ratio), documented fixed-tau1 scaling', note='Additional gain cases are read from EXTRA_stuff; q companions are not needed for dipole pulse metrics. Resonant frequency headers and normalized drive reconstructed from dipole/alpha are checked. Run settings inherited from original figure; output repeatability is not a timestep convergence test. Extra/1.12Gthq.dat appears mislabeled and is not used.', files=rows),indent=2))
    subprocess.run(['gnuplot',str(BASE/'gnuplot'/ (STEM+'.gp'))],check=True)
    print(f'Rendered {len(rows)} gain cases')

if __name__ == '__main__':
    main()
