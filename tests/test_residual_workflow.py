#!/usr/bin/env python3
"""Integration regression; pass the compiled emission_convergence driver."""
import argparse
import math
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]

def comments(path):
    return dict(line[2:].split(' = ',1) for line in path.read_text().splitlines() if line.startswith('# ') and ' = ' in line)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver',type=Path)
    parser.add_argument('--reference-driver',type=Path,help='optional pre-migration driver for prepare-mode metadata comparison')
    args=parser.parse_args();binary=args.driver.resolve()
    source=ROOT/'data/output/convergence_validation/inputs'
    with tempfile.TemporaryDirectory(prefix='residual-workflow-') as temp:
        base=Path(temp)
        def run(name,mode,factor,*extra):
            folder=base/name;folder.mkdir()
            with (folder/'run.log').open('w') as log:
                subprocess.run([str(binary),str(source/'input08.dat'),str(source/'time.dat'),str(folder),mode,str(factor),*map(str,extra)],check=True,stdout=log,stderr=log,timeout=120)
            return folder
        saturation=[]
        for gain in (1.05,1.10,1.15):
            folder=run(f'normalization-{gain}','prepare',1,'threshold',gain,1e-9)
            md=dict(line.split('\t',1) for line in (folder/'driver_metadata.tsv').read_text().splitlines())
            sat=float(md['Esat_solver_units']);plasma=float(md['plasma_energy_eV'])
            assert md['saturation_scale_convention']=='inversion_one_plasma_scaled_saturation_v1'
            assert md['tau1_policy']=='require_match'
            assert md['supplied_tau1_matches']=='true' and md['tau1_mismatch_overridden']=='false'
            assert float(md['tau1_over_tau2'])==float(md['effective_tau1_over_tau2'])==5.
            tau1=10.*plasma/float(md['linewidth_eV'])
            for key in ('tau1_normalized','energy_derived_tau1_normalized','supplied_tau1_normalized'):
                assert math.isclose(float(md[key]),tau1,rel_tol=1e-14),key
            assert math.isclose(float(md['tau2_normalized']),tau1/5.,rel_tol=1e-14)
            assert math.isclose(float(md['Isat_solver_units']),1./(abs(float(md['gain']))*tau1),rel_tol=1e-14)
            if args.reference_driver:
                reference=base/f'reference-{gain}';reference.mkdir()
                subprocess.run([str(args.reference_driver.resolve()),str(source/'input08.dat'),str(source/'time.dat'),str(reference),'prepare','1','threshold',str(gain),'1e-9'],check=True,stdout=subprocess.DEVNULL,timeout=30)
                old=dict(line.split('\t',1) for line in (reference/'driver_metadata.tsv').read_text().splitlines())
                for key,value in old.items():
                    try:
                        expected=float(value)
                    except ValueError:
                        assert md[key]==value,key
                    else:
                        assert math.isclose(float(md[key]),expected,rel_tol=1e-14,abs_tol=0.),key
            assert math.isclose(float(md['E0_over_Esat']),1e-9,rel_tol=1e-14)
            assert math.isclose(float(md['raw_E0_V_m'])/float(md['Esat_V_m']),1e-9,rel_tol=1e-14)
            assert math.isclose(sat*float(md['field_unit_V_m']),float(md['Esat_V_m']),rel_tol=1e-12)
            assert math.isclose(sat*math.sqrt(plasma),float(md['legacy_Esat_solver_units']),rel_tol=1e-12)
            saturation.append(float(md['Esat_V_m']))
        assert max(saturation)/min(saturation)-1<1e-12
        baseline=None
        for name,mode,factor in [('fixed','fixed',1),('automatic','automatic',1),('tight','automatic',.5)]:
            folder=run(name,mode,factor)
            expected=ROOT/'data/output/convergence_validation'/name/'08'/'numtime.dat'
            assert (folder/'numtime.dat').read_bytes()==expected.read_bytes(),name
            if name=='automatic': baseline=folder
        explicit_zero=run('explicit-zero','automatic',1,2.811,'input','input',.1,0.0)
        assert (explicit_zero/'numtime.dat').read_bytes()==(baseline/'numtime.dat').read_bytes()
        observed=[]
        for cadence in (.1,.2):
            folder=run(f'observation-{cadence}','automatic',1,2.811,'input','input',cadence,40,3,.05,1e-8,1e-4)
            md=comments(folder/'numtime.dat');original=comments(baseline/'numtime.dat')
            assert all(md[k]==v for k,v in original.items() if k.startswith('convergence_'))
            assert md['post_acceptance_status']=='complete'
            text=(folder/'numtime.dat').read_text()
            if cadence==.1:
                old=[x for x in (baseline/'numtime.dat').read_text().splitlines() if x and not x.startswith('#')]
                new=[x for x in text.splitlines() if x and not x.startswith('#')]
                assert old[:-1]==new[:len(old)-1], 'trajectory prefix changed'
                assert float(new[-1].split()[0])>float(old[-1].split()[0])+120
            observed.append((folder/'diagnostic_windows.tsv').read_bytes())
        assert observed[0]==observed[1], 'diagnostic depends on output cadence'
    print('SI normalization at all three gains, legacy disabled behavior, same-trajectory prefix and diagnostic cadence independence: PASS')

if __name__=='__main__': main()
