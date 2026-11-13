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
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--timeout-seconds', type=float, default=120)
    p.add_argument('--cadence-ps', type=float, default=.1)
    a = p.parse_args()
    for value in (a.timeout_seconds, a.cadence_ps):
        if not math.isfinite(value) or value <= 0:
            p.error('timeout and cadence must be finite and positive')
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
    sources = [ROOT/'src/emission_convergence.cxx', Path(__file__).resolve()]
    manifest = dict(frequencies=frequencies, input_sha256=digest(a.input), time_sha256=digest(a.time),
                    repo_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
                    submodule_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ngm, text=True).strip(),
                    source_sha256={str(x.relative_to(ROOT)):digest(x) for x in sources},
                    cadence_ps=a.cadence_ps, timeout_seconds=a.timeout_seconds,
                    normalization='src/Esat.cxx: Esat=sqrt(Dome/(10*abs(G))); manuscript agreement not independently verified',
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
                command=[str(binary),str(a.output/'input.dat'),str(a.output/'time.dat'),str(folder),'automatic',str(factor),str(frequency),gain,'1e-9',str(a.cadence_ps)]
                (folder/'command.json').write_text(json.dumps(command)+'\n')
                for stale in ('numtime.dat', 'driver_metadata.tsv'):
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
                md=metadata(folder/'numtime.dat'); row.update(md)
                reason=md.get('termination_reason') if status=='completed' else None
                row['classification']={'sustained_flat':'flat_tolerance_dependent','sustained_periodic':'periodic'}.get(reason,'unclassified')
                row['relative_peak_to_peak_modulation']=None
                if row['classification']!='unclassified':
                    mean=float(md['convergence_mean_intensity']); pp=float(md['convergence_peak_to_peak'])
                    row['relative_peak_to_peak_modulation']=pp/mean if mean>0 else None
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
        writer.writeheader(); writer.writerows(rows)

if __name__=='__main__':
    main()
