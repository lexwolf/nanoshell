"""Numerical decomposition and report-support diagnostics; no nonlinear runs."""
from pathlib import Path
import sys,csv,json
import numpy as np
from scipy.linalg import eig
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
mod=__import__('[OMEx]_turn_on_linear');ROOT=mod.ROOT;OUT=mod.OUT;write=mod.write;models=mod.models();eth=2.8121971959422813
IMG=ROOT/'img/output/[OMEx]_turn_on_mechanisms';OLD=ROOT/'data/output/gain_switching_corrected_20261114'
read=lambda n:list(csv.DictReader((OUT/'tables'/n).open()))
linear={float(r['energy_eV']):{k:float(v) for k,v in r.items()} for r in read('linear_summary.csv')};tra=read('trajectory_comparison.csv');forcing=[]
for e,m in models.items():
 A=m['A0']+m['D'];lam,V=eig(A);W=np.linalg.inv(V);j=int(lam.real.argmax());v=V[:,j];w=W[j];l=lam[j];b=m['b0']+m['db'];cv=m['c'][:3]@v
 forcing.append(dict(energy_eV=e,abs_lambda_ps_inv=abs(l),modal_forcing_output_rate_over_Esat=abs(cv*(w@b))/m['Esat'],growth_sensitivity_to_N_ps_inv=float(np.real(w@m['D']@v)),frozen_seed_over_Esat=abs(cv*(w@b)/l)/m['Esat'],ramp_seed_over_Esat=linear[e]['effective_output_seed_over_Esat'],time_to_instability_undepleted_ps=-m['tau1']*np.log1p(-linear[e]['instability_N'])))
write('forcing_projection.csv',forcing)
L=next(r for r in tra if r['case']=='pair_left_dt1');R=next(r for r in tra if r['case']=='pair_right_dt1')
sL=float(L['sigma_ps_inv']);sR=float(R['sigma_ps_inv']);uL=np.log(float(L['effective_seed_over_Esat']));uR=np.log(float(R['effective_seed_over_Esat']));vL=np.log(np.sqrt(.81*float(L['maximum_I'])));vR=np.log(np.sqrt(.81*float(R['maximum_I'])))
rate=((vL+vR-uL-uR)/2)*(1/sR-1/sL);seed=-(uR-uL)*(1/sR+1/sL)/2;criterion=(vR-vL)*(1/sR+1/sL)/2
report=dict(observed_delta_ps=float(R['t_on_own81_ps'])-float(L['t_on_own81_ps']),linear_delta_ps=float(R['t_on_own81_linear_prediction_ps'])-float(L['t_on_own81_linear_prediction_ps']),symmetric_growth_rate_contribution_ps=rate,symmetric_seed_contribution_ps=seed,symmetric_own_peak_level_contribution_ps=criterion,remaining_nonlinear_approach_contribution_ps=float(R['nonlinear_criterion_correction_ps'])-float(L['nonlinear_criterion_correction_ps']),note='Algebraic decomposition of asymptotic linear extrapolation; remaining term measured from nonlinear crossing, not a separately simulated mechanism.')
(OUT/'tables/pair_decomposition.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
d=1e-6;sigma=linear[eth]['sigma_ps_inv'];nup=(linear[eth+d]['imag_ps_inv']-linear[eth-d]['imag_ps_inv'])/(2*d)
local={'amplitude_sigma_ps_inv':sigma,'d_nu_d_energy_ps_inv_per_eV':nup,'modal_forcing_width_eV':sigma/abs(nup),'predicted_even_quadratic_coefficient_ps_per_eV2':nup*nup/(2*sigma**3),'assumptions':'Locally freeze sigma, modal forcing numerator, output criterion and inversion-ramp correction; vary imaginary pole detuning. Diagnostic approximation, not an exact nonlinear law.'}
(OUT/'tables/central_pole_prediction.json').write_text(json.dumps(local,indent=2))
# Smoothness at the smallest detunings, checked at BOTH integration steps.
centr=[]
for prefix in ['central_','half_']:
 c=next(r for r in tra if r['case']==prefix+'+0e+00');t0=float(c['t_on_own81_ps'])
 for d in [1e-3,1e-5,1e-6]:
  a=next(r for r in tra if r['case']==prefix+f'{d:+.0e}');b=next(r for r in tra if r['case']==prefix+f'{-d:+.0e}')
  centr.append(dict(step='baseline' if prefix=='central_' else 'half',abs_detuning_eV=d,even_delay_ps=(float(a['t_on_own81_ps'])+float(b['t_on_own81_ps']))/2-t0,odd_delay_ps=(float(a['t_on_own81_ps'])-float(b['t_on_own81_ps']))/2))
write('central_step_smoothness.csv',centr)
saturation=[]
for name in ['pair_left_dt1','pair_right_dt1']:
 row=next(r for r in tra if r['case']==name);e=float(row['energy_eV']);m=models[e];Ncrit=linear[e]['instability_N'];vals=eig(m['A0']+Ncrit*m['D'],right=False);v=vals[np.argmin(abs(vals.real))];tr=np.loadtxt(ROOT/row['source']/'trajectory.tsv');tail=tr[tr[:,0]>=tr[-1,0]-10]
 saturation.append(dict(case=name,instability_N=Ncrit,tail_mean_N=float(tail[:,4].mean()),tail_mean_I=float(tail[:,1].mean()),critical_mode_sigma_ps_inv=v.real,critical_mode_nu_ps_inv=v.imag,critical_mode_lab_energy_eV=e-(9.6/m['k'])*v.imag))
write('pair_saturation.csv',saturation)
# Denser N discretization cross-check of modal continuation has no crossings in this set.
checks=[]
for e in [eth,2.5309774763480535,3.1201269537068566]:
 m=models[e];from itertools import permutations
 for count in [501,1001]:
  prev=None;minsep=np.inf
  for N in np.linspace(0,1,count):
   val=eig(m['A0']+N*m['D'],right=False);order=np.argsort(val.imag) if prev is None else min(permutations(range(3)),key=lambda p:sum(abs(val[list(p)]-prev)**2));val=val[list(order)];prev=val;minsep=min(minsep,min(abs(val[i]-val[j]) for i in range(3) for j in range(i)))
  checks.append(dict(energy_eV=e,N_samples=count,dominant_mode=int(val.real.argmax()),sigma_ps_inv=max(val.real),min_separation_ps_inv=minsep))
write('mode_tracking_checks.csv',checks)
# Forced-propagation check at exactly aligned observation times, avoiding intensity-zero divisions.
checks=[]
for name in ['pair_left_dt1','pair_right_dt1']:
 row=next(r for r in tra if r['case']==name);m=models[float(row['energy_eV'])];folder=ROOT/row['source'];tr=np.loadtxt(folder/'trajectory.tsv');tr=tr[tr[:,0]<=4];tt=tr[:,0]+float(row['step_ps'])
 def rhs(t,q):
  N=-np.expm1(-t/m['tau1']);return (m['A0']+N*m['D'])@q+m['b0']+N*m['db']
 sol=solve_ivp(rhs,[0,tt[-1]],np.zeros(3,complex),t_eval=tt,method='DOP853',rtol=2e-11,atol=1e-25);assert sol.success
 dp=m['c'][:3]@sol.y+m['c'][3]*m['E0'];meas=tr[:,2]+1j*tr[:,3]
 checks.append(dict(case=name,time_interval_ps=[0,4],relative_L2_complex_dipole_error=float(np.linalg.norm(meas-dp)/np.linalg.norm(dp)),comparison='continuous ODE linear, undepleted inversion, vs split-step nonlinear; native historical labels shifted forward by dt'))
(OUT/'tables/forced_validation.json').write_text(json.dumps(checks,indent=2))
plt.rcParams.update({'font.size':10,'figure.constrained_layout.use':True})
fig,ax=plt.subplots(2,2,figsize=(10,7))
for name,label in [('central_+0e+00','Exact threshold'),('central_+1e-02','Threshold +10 meV'),('residual_default_g1.1_0050','Left interior')]:
 found=[r for r in tra if r['case']==name]
 if not found:continue
 row=found[0];tr=np.loadtxt(ROOT/row['source']/'trajectory.tsv');e=float(row['energy_eV']);z=np.load(OUT/'tables'/f'linear_{e:.17g}.npz')
 ax[0,0].semilogy(tr[:,0],tr[:,1],label=label);ax[0,1].plot(tr[:,0],tr[:,4],label=label)
 if name=='central_+0e+00':
  ax[1,0].semilogy(z['t'],z['I'],label='Inversion starts at zero');ax[1,0].semilogy(z['t'],z['forced_N1_I'],'--',label='Inversion fixed at one')
for a in ax.flat:a.legend(fontsize=8) if a.get_legend_handles_labels()[0] else None
ax[0,0].set(xlim=(0,15),ylim=(1e-17,.1),xlabel='Time (ps)',ylabel='I / Isat');ax[0,1].set(xlim=(0,15),xlabel='Time (ps)',ylabel='Inversion N');ax[1,0].set(xlim=(0,.4),ylim=(1e-19,1e-12),xlabel='Time (ps)',ylabel='Forced linear I / Isat')
f=sorted(forcing,key=lambda r:r['energy_eV']);ax[1,1].semilogy([r['energy_eV'] for r in f],[r['ramp_seed_over_Esat'] for r in f],label='Ramped-inversion seed');ax[1,1].semilogy([r['energy_eV'] for r in f],[r['frozen_seed_over_Esat'] for r in f],'--',label='N=1 seed');ax[1,1].legend(fontsize=8);ax[1,1].set(xlabel='Probe energy (eV)',ylabel='Effective dipolar seed / Esat')
for ext in ['pdf','png']:fig.savefig(IMG/f'[OMEx]_growth_and_inversion.{ext}',dpi=180)
plt.close(fig)
