"""Gate nonlinear interpretation on passive response, carrier covariance and step refinement."""
from pathlib import Path
import json,csv
import numpy as np
m=__import__('[OMEx]_metal_envelope_validation');p=m.OUT;c=m.constants();eth=c['old_threshold_eV'];ep=c['plasma_eV'];rows=[]
for tag in ['left','center','right']:
 f=p/'runs'/('passive_'+tag);r=json.loads((f/'metrics.json').read_text());a=np.loadtxt(f/'trajectory.tsv');a=a[a[:,0]>1.5];d=a[:,2]+1j*a[:,3];response=d/r['E0']*np.exp(-1j*(r['reference_eV']-r['probe_eV'])/ep*c['k_ps_inv']*a[:,0]);e=r['probe_eV'];w=e/ep;epsm=c['eps_inf']-1/(w*(w+2j*c['gamma_normalized']));eb=c['eps_core'];ee=c['eps_host'];r3=.6**3;analytic=((epsm-ee)*(eb+2*epsm)+r3*(eb-epsm)*(ee+2*epsm))/((epsm+2*ee)*(eb+2*epsm)+2*r3*(epsm-ee)*(eb-epsm));err=abs(response.mean()-analytic)/abs(analytic);assert err<1e-7;rows.append(dict(case=tag,probe_eV=e,relative_complex_passive_response_error=err))
m.write('passive_time_validation.csv',rows)
cov=[]
for dt in [.125,.0625]:
 a=np.loadtxt(p/'runs'/f'threshold_dt{dt:g}'/'trajectory.tsv');b=np.loadtxt(p/'runs'/f'carrier_dt{dt:g}'/'trajectory.tsv');assert np.array_equal(a[:,0],b[:,0]);da=a[:,2]+1j*a[:,3];db=(b[:,2]+1j*b[:,3])*np.exp(-1j*.1/ep*c['k_ps_inv']*b[:,0]);err=np.linalg.norm(da-db)/np.linalg.norm(da);relI=np.max(abs(a[:,1]-b[:,1]))/max(a[:,1]);Nerr=float(np.max(abs(a[:,4]-b[:,4])));assert err<1e-7 and relI<1e-7 and Nerr<1e-8;cov.append(dict(dt=dt,reference_shift_eV=.1,relative_complex_dipole_L2_error=err,intensity_max_error_over_peak=relI,inversion_max_absolute_error=Nerr))
m.write('carrier_time_validation.csv',cov)
steps=[]
for tag in ['threshold','left','right']:
 fine=json.loads((p/'runs'/f'{tag}_dt0.0625'/'metrics.json').read_text())
 for dt in [1,.5,.25,.125,.0625]:
  f=p/'runs'/f'{tag}_dt{dt:g}'/'metrics.json'
  if not f.exists():continue
  a=json.loads(f.read_text());steps.append(dict(case=tag,dt_normalized=dt,ton81_ps=a['ton81_ps'],difference_from_finest_ps=a['ton81_ps']-fine['ton81_ps'],peak_relative_difference=a['maximum_J']/fine['maximum_J']-1))
m.write('step_convergence.csv',steps)
assert max(abs(x['difference_from_finest_ps']) for x in steps if x['dt_normalized']<=.25)<1e-6
summary={'passed':True,'primary_dt_normalized':.125,'finest_dt_normalized':.0625,'max_primary_timing_error_ps':max(abs(x['difference_from_finest_ps']) for x in steps if x['dt_normalized']==.125),'max_passive_time_relative_error':max(x['relative_complex_passive_response_error'] for x in rows),'max_carrier_time_L2_error':max(x['relative_complex_dipole_L2_error'] for x in cov)}
(p/'tables/pilot_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
