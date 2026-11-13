#!/usr/bin/env python3
"""Isolated residual scan; all convergence decisions/statistics come from CUP."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
GAINS = ('1.05', '1.10', '1.15')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def metadata(path):
    result = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith('# ') and ' = ' in line:
                key, value = line[2:].split(' = ', 1)
                result[key] = value
    return result

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frequencies', required=True, type=Path, help='JSON: each of 1.05, 1.10, 1.15 maps to an explicit list of eV values or "threshold"; empty lists allowed')
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--time', required=True, type=Path)
    p.add_argument('--manuscript', type=Path, help='optional manuscript source snapshot for provenance')
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--timeout-seconds', type=float, default=120)
    p.add_argument('--cadence-ps', type=float, default=.1)
    p.add_argument('--diagnostic-window-ps', type=float, default=40., help='0 disables continuation; otherwise minimum duration per late window')
    p.add_argument('--diagnostic-windows', type=int, default=3)
    p.add_argument('--range-rtol', type=float, default=.05)
    p.add_argument('--range-atol', type=float, default=1e-8)
    p.add_argument('--mean-rtol', type=float, default=1e-4)
    a = p.parse_args()
    for value in (a.timeout_seconds, a.cadence_ps):
        if not math.isfinite(value) or value <= 0:
            p.error('timeout and cadence must be finite and positive')
    if not math.isfinite(a.diagnostic_window_ps) or not 0 <= a.diagnostic_window_ps <= 100000 or not 3 <= a.diagnostic_windows <= 4096:
        p.error('diagnostic duration must be in [0,100000] ps; count must be in [3,4096]')
    if any(not math.isfinite(x) or x <= 0 for x in (a.range_rtol, a.mean_rtol)) or not math.isfinite(a.range_atol) or a.range_atol < 0:
        p.error('invalid diagnostic comparison tolerances')
    frequencies = json.loads(a.frequencies.read_text())
    if set(frequencies) != set(GAINS):
        p.error('provide all three gain keys (empty lists allowed)')
    for values in frequencies.values():
        if not isinstance(values, list) or len(set(map(str, values))) != len(values):
            p.error('frequency lists must be lists without duplicates')
        for value in values:
            if value != 'threshold' and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0):
                p.error('frequencies must be positive finite numbers in eV or "threshold"')
    a.output = a.output.resolve()
    a.output.mkdir(parents=True, exist_ok=True)
    binary = a.output/'emission_convergence'
    ngm = ROOT/'extern/nano_geo_matrix'
    sources = [ROOT/'src/emission_convergence.cxx', Path(__file__).resolve(), ROOT/'src/residual_diagnostic.hpp', ROOT/'docs/residual_modulation.md']
    manifest = dict(frequencies=frequencies, input_sha256=digest(a.input), time_sha256=digest(a.time),
                    repo_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
                    submodule_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ngm, text=True).strip(),
                    source_sha256={str(x.relative_to(ROOT)):digest(x) for x in sources},
                    cadence_ps=a.cadence_ps, timeout_seconds=a.timeout_seconds,
                    normalization='manuscript_SI_to_plasma_scaled_CUP_v1: Esat=sqrt(Dome/(10*abs(G)*Ome_p))',
                    manuscript_sha256=digest(a.manuscript) if a.manuscript else None,
                    diagnostic=dict(window_ps=a.diagnostic_window_ps, windows=a.diagnostic_windows, range_rtol=a.range_rtol, range_atol=a.range_atol, mean_rtol=a.mean_rtol),
                    repo_diff_sha256=hashlib.sha256(subprocess.check_output(['git','diff','HEAD'], cwd=ROOT)).hexdigest(),
                    submodule_diff_sha256=hashlib.sha256(subprocess.check_output(['git','diff','HEAD'], cwd=ngm)).hexdigest(),
                    build='g++ -std=c++17 -O2; gsl gslcblas armadillo')
    manifest_path=a.output/'manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        p.error('resume configuration/provenance mismatch; use a new output directory')
    (a.output/'repository.patch').write_bytes(subprocess.check_output(['git','diff','HEAD'], cwd=ROOT))
    (a.output/'submodule.patch').write_bytes(subprocess.check_output(['git','diff','HEAD'], cwd=ngm))
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    (a.output/'input.dat').write_bytes(a.input.read_bytes())
    (a.output/'time.dat').write_bytes(a.time.read_bytes())
    (a.output/'frequencies.json').write_bytes(a.frequencies.read_bytes())
    if a.manuscript:
        (a.output/'manuscript.tex').write_bytes(a.manuscript.read_bytes())
    for source in sources:
        (a.output/source.name).write_bytes(source.read_bytes())
    with (a.output/'build.log').open('w') as log:
        subprocess.run(['g++','-std=c++17','-O2',*[f'-I{ngm/x}' for x in ('include','modules','modules/cup')], '-I/usr/include/eigen3',str(sources[0]),'-o',str(binary),'-lgsl','-lgslcblas','-lm','-larmadillo'], stdout=log, stderr=log, check=True)
    rows=[]
    for gain in GAINS:
        for index, frequency in enumerate(frequencies[gain]):
            for factor in (1., .5):
                folder=a.output/f'gain_{gain}'/f'frequency_{index:04d}'/f'tolerance_{factor:g}'
                folder.mkdir(parents=True, exist_ok=True)
                record=folder/'result.json'
                if record.exists() and json.loads(record.read_text())['run_status']=='completed':
                    rows.append(json.loads(record.read_text())); continue
                command=[str(binary),str(a.output/'input.dat'),str(a.output/'time.dat'),str(folder),'automatic',str(factor),str(frequency),gain,'1e-9',str(a.cadence_ps),str(a.diagnostic_window_ps),str(a.diagnostic_windows),str(a.range_rtol),str(a.range_atol),str(a.mean_rtol)]
                (folder/'command.json').write_text(json.dumps(command)+'\n')
                for stale in ('numtime.dat', 'driver_metadata.tsv', 'diagnostic_metadata.tsv', 'diagnostic_windows.tsv', 'accepted_metadata.dat'):
                    (folder/stale).unlink(missing_ok=True)
                with (folder/'run.log').open('w') as log:
                    try:
                        process=subprocess.run(command, stdout=log, stderr=log, timeout=a.timeout_seconds)
                        status='completed' if process.returncode==0 else 'failed'
                    except subprocess.TimeoutExpired:
                        status='wall_timeout'
                row=dict(gain_ratio=gain, requested_frequency_eV=frequency, tolerance_factor=factor, run_status=status, folder=str(folder.relative_to(a.output)))
                driver=folder/'driver_metadata.tsv'
                if driver.exists():
                    row.update(dict(line.split('\t',1) for line in driver.read_text().splitlines()))
                md=metadata(folder/'accepted_metadata.dat'); md.update(metadata(folder/'numtime.dat')); row.update(md)
                reason=md.get('accepted_termination_reason') or (md.get('termination_reason') if status=='completed' else None)
                row['classification']={'sustained_flat':'flat_tolerance_dependent','sustained_periodic':'periodic'}.get(reason,'unclassified')
                row['accepted_relative_interval_range']=None
                row['late_variation_assessment']='unresolved'
                row['diagnostic_windows']=[]
                if row['classification']!='unclassified':
                    mean=float(md['convergence_mean_intensity']); pp=float(md['convergence_peak_to_peak'])
                    row['accepted_relative_interval_range']=pp/mean if mean>0 else None
                    row['accepted_interval_peak_to_peak']=pp
                    row['accepted_mean_I_over_Isat']=mean/float(row['Isat_solver_units'])
                diag=folder/'diagnostic_metadata.tsv'
                if diag.exists():
                    row.update(dict(line.split('\t',1) for line in diag.read_text().splitlines()))
                    row['late_variation_assessment']=row['diagnostic_status'] if status=='completed' else 'unresolved'
                    if status!='completed': row['diagnostic_complete']='0'
                    with (folder/'diagnostic_windows.tsv').open() as stream:
                        for window in csv.DictReader(stream, delimiter='\t'):
                            row['diagnostic_windows'].append({k:(None if v=='null' else float(v)) for k,v in window.items()})
                if 'actual_output_interval_ps' in md:
                    row['integration_step_ps']=float(md['actual_output_interval_ps'])/int(md['output_interval_steps'])
                record.write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
                rows.append(row)
                # Preserve a usable partial summary after every point.
                save(a.output, rows)
    save(a.output, rows)

def save(output, rows):
    (output/'summary.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with (output/'summary.tsv').open('w') as stream:
        writer=csv.DictWriter(stream, fieldnames=keys, delimiter='\t')
        writer.writeheader()
        writer.writerows({k:(json.dumps(v) if isinstance(v,list) else v) for k,v in r.items()} for r in rows)
    with (output/'late_windows.tsv').open('w') as stream:
        fields=['gain_ratio','requested_frequency_eV','tolerance_factor','start_ps','end_ps','mean_intensity','minimum_intensity','maximum_intensity','peak_to_peak','relative_interval_range','cycles']
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t'); writer.writeheader()
        for row in rows:
            for window in row.get('diagnostic_windows',[]):
                writer.writerow({**{k:row[k] for k in fields[:3]},**window})

if __name__=='__main__':
    main()
