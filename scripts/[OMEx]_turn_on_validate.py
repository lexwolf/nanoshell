"""Focused checks of modal normalization, propagation, refinement and reuse."""
from pathlib import Path
import ast,csv,json,subprocess
import numpy as np
from scipy.linalg import eig,expm
mod=__import__('[OMEx]_turn_on_linear');OUT=mod.OUT;ROOT=mod.ROOT
read=lambda n:list(csv.DictReader((OUT/'tables'/n).open()))
checks={}
a=np.genfromtxt(OUT/'tables/tracked_modes.csv',names=True,delimiter=',');checks['minimum_mode_separation_ps_inv']=float(a['min_eigenvalue_separation_ps_inv'].min());checks['max_eigenvector_condition']=float(a['eigenvector_matrix_condition'].max());assert checks['minimum_mode_separation_ps_inv']>90
checks['max_biorthogonal_residual']=max(float(r['biorthogonal_residual']) for r in read('linear_summary.csv'));assert checks['max_biorthogonal_residual']<1e-12
errors=[]
for e,m in mod.models().items():
 A=m['A0']+m['D'];b=m['b0']+m['db'];v,V=eig(A);W=np.linalg.inv(V)
 for t in [.01,.1,1.]:
  direct=np.linalg.solve(A,(expm(A*t)-np.eye(3))@b);modal=V@((np.exp(v*t)-1)/v*(W@b));errors.append(np.linalg.norm(direct-modal)/np.linalg.norm(direct))
checks['max_frozen_forced_modal_vs_expm_relative_error']=float(max(errors));assert max(errors)<1e-9
r=read('trajectory_comparison.csv');checks['max_relative_fitted_sigma_error']=max(abs(float(x['field_sigma_fit_ps_inv'])/float(x['sigma_ps_inv'])-1) for x in r);assert checks['max_relative_fitted_sigma_error']<.001
checks['max_fixed_threshold_prediction_error_ps']=max(abs(float(x['fixed_I_1e_4_crossing_ps'])-float(x['fixed_I_1e_4_linear_prediction_ps'])) for x in r);assert checks['max_fixed_threshold_prediction_error_ps']<.01
smooth=read('central_step_smoothness.csv');checks['max_half_step_even_delay_relative_change']=max(abs(float(next(x for x in smooth if x['step']=='half' and x['abs_detuning_eV']==a['abs_detuning_eV'])['even_delay_ps'])/float(a['even_delay_ps'])-1) for a in smooth if a['step']=='baseline');assert checks['max_half_step_even_delay_relative_change']<.001
for case in ['central_+0e+00','central_+1e-03','central_-1e-03']:
 rr=[float(x['own81_crossing_ps']) for x in read('horizon_checks.csv') if x['case']==case and float(x['horizon_ps'])>=20];assert np.ptp(rr)<1e-9
rounded=next(x for x in r if x['case']=='central_rounded');old=next(x for x in r if x['case']=='residual_default_g1.1_0099');checks['rounded_reproduction_delay_difference_ps']=float(rounded['t_on_own81_ps'])-float(old['t_on_own81_ps']);assert abs(checks['rounded_reproduction_delay_difference_ps'])<1e-6
forced=json.loads((OUT/'tables/forced_validation.json').read_text());checks['max_forced_continuous_vs_split_step_complex_dipole_L2_error']=max(x['relative_L2_complex_dipole_error'] for x in forced);assert checks['max_forced_continuous_vs_split_step_complex_dipole_L2_error']<.001
for p in (ROOT/'scripts').glob('[[]OMEx]_turn_on_*.py'):ast.parse(p.read_text())
subprocess.run(['bash','-n',str(ROOT/'scripts/[OMEx]_turn_on_mechanisms.bash')],check=True)
checks['all_checks_passed']=True;(OUT/'tables/validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
