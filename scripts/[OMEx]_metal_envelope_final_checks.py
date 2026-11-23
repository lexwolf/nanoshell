"""Final targeted checks: independent frequency-domain threshold and observed controls."""
import json,csv
import numpy as np
from scipy.optimize import root
m=__import__('[OMEx]_metal_envelope_validation');p=m.OUT;c=m.constants();s=json.loads((p/'tables/thresholds.json').read_text());ep=c['plasma_eV'];wg=c['old_threshold_eV']/ep
# Independent constitutive/boundary denominator, rather than eigenvalue zero.
def equations(x):
 G,e=x;w=e/ep;em=c['eps_inf']-1/(w*(w+2j*c['gamma_normalized']));ec=c['eps_core']+G/(c['tau2']*(w-wg)+1j);ee=c['eps_host'];D=(em+2*ee)*(ec+2*em)+2*.6**3*(ec-em)*(em-ee);return [D.real,D.imag]
r=root(equations,[s['fixed_gain_center_Gth'],s['fixed_gain_center_threshold_emission_eV']],tol=1e-11);assert r.success
assert abs(r.x[0]-s['fixed_gain_center_Gth'])<1e-10 and abs(r.x[1]-s['fixed_gain_center_threshold_emission_eV'])<1e-9
rows=list(csv.DictReader((p/'tables/trajectory_comparison.csv').open()));cur={x['case']:x for x in rows}
for q in (p/'runs').glob('*'):
 if (q/'result.json').exists():assert json.loads((q/'result.json').read_text())['returncode']==0
# Common absolute-intensity crossing preserves the near-center minimum and the original-pair small asymmetry.
assert float(cur['threshold_dt0.125']['fixed_J_1e_4_ps'])<float(cur['central_+1e-03']['fixed_J_1e_4_ps'])
assert float(cur['threshold_dt0.125']['fixed_J_1e_4_ps'])<float(cur['central_-1e-03']['fixed_J_1e_4_ps'])
checks=[]
for d in [-.001,.001]:
 a=cur[f'central_{d:+.0e}'];b=cur[f'centralhalf_{d:+.0e}'];checks.append(dict(detuning_eV=d,half_minus_primary_ton_ps=float(b['ton81_ps'])-float(a['ton81_ps'])));assert abs(checks[-1]['half_minus_primary_ton_ps'])<1e-7
m.write('central_step_checks.csv',checks)
# Evaluate horizon changes rather than assuming trajectory maximum is final-state intensity.
h=list(csv.DictReader((p/'tables/horizon_checks.csv').open()));horizon=[]
for name in ['threshold_dt0.125','left_dt0.125','right_dt0.125']:
 a=[float(x['ton81_ps']) for x in h if x['case']==name];span=np.ptp(a);horizon.append(dict(case=name,ton_span_20_to_30ps=span));assert span<1e-5
m.write('primary_horizon_check.csv',horizon)
validation=dict(passed=True,independent_constitutive_Gth=float(r.x[0]),independent_constitutive_threshold_eV=float(r.x[1]),eigenvalue_Gth=s['fixed_gain_center_Gth'],threshold_residual=equations(r.x),run_count=sum((q/'result.json').exists() for q in (p/'runs').glob('*')),max_central_half_step_timing_difference_ps=max(abs(x['half_minus_primary_ton_ps']) for x in checks),common_absolute_crossing_minimum_preserved=True)
(p/'tables/final_validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2))
