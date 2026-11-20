"""Compare tracked modes, forced linear evolution and native-step nonlinear metrics."""
from pathlib import Path
import csv,json,sys,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/turn_on_mechanisms_20261117';OLD=ROOT/'data/output/gain_switching_corrected_20261114';IMG=ROOT/'img/output/[OMEx]_turn_on_mechanisms';IMG.mkdir(exist_ok=True)
native=__import__('[OMEx]_native_io').load
linear=__import__('[OMEx]_turn_on_linear');write=linear.write;models=linear.models()
eth=2.8121971959422813
read=lambda p:dict(x.split('\t',1) for x in p.read_text().splitlines())
ls=list(csv.DictReader((OUT/'tables/linear_summary.csv').open()));lin={float(r['energy_eV']):{k:float(v) for k,v in r.items()} for r in ls}
selection=json.loads((OUT/'inputs/reused_selection.json').read_text());cases=[]
for x in selection:cases.append((x['run'],OLD/'runs'/x['run']))
for name in ['pair_left_dt1','pair_right_dt1','pair_left_dt0.5','pair_right_dt0.5']:cases.append((name,OLD/'runs'/name))
cases += [(q.name,q) for q in sorted((OUT/'runs').glob('*')) if (q/'metrics.tsv').exists()]
rows=[];horizons=[];series={};central=[];inputs=[]
def cross(I,dt,level):
 hit=np.flatnonzero((I[:-1]<level)&(I[1:]>=level))
 if len(hit)==0:return np.nan
 i=hit[0];return (i+(level-I[i])/(I[i+1]-I[i]))*dt
for name,folder in cases:
 par=read(folder/'parameters.tsv');met=read(folder/'metrics.tsv');e=float(par['frequency_eV']);dt=float(par['integration_step_ps']);l=lin[e];sig=l['sigma_ps_inv'];seed=l['effective_output_seed_over_Esat']
 assert abs(float(par['gain_ratio'])-1.1)<1e-12 and par['pump']=='step' and float(par['initial_N'])==0 and par['initial_q']=='0,0,0'
 assert abs(float(par['E0_over_Esat'])/1e-9-1)<1e-12 and abs(float(par['Esat'])/models[e]['Esat']-1)<1e-12
 I=native(folder);t=np.arange(len(I))*dt;tr=np.loadtxt(folder/'trajectory.tsv');series[name]=(tr,I[::max(1,int(.002/dt))].copy(),dt*max(1,int(.002/dt)))
 peak=float(I.max());ton=cross(I,dt,.81*peak)
 # Native intensities, after inversion ramp and well before appreciable depletion.
 lo=1e-6;hi=min(1e-5,.001*peak);mask=(I>=lo)&(I<=hi)&(t>=max(1.,3/sig))
 if mask.sum()>=20:
  x=t[mask];y=np.log(I[mask]);coef=np.polyfit(x,y,1);res=y-np.polyval(coef,x);r2=1-np.sum(res**2)/np.sum((y-y.mean())**2)
  fit=dict(fit_t_start_ps=float(x.min()),fit_t_end_ps=float(x.max()),fit_samples=int(mask.sum()),intensity_log_slope_ps_inv=float(coef[0]),field_sigma_fit_ps_inv=float(coef[0]/2),fit_R2=float(r2),fit_inversion_min=float(np.interp(x,tr[:,0],tr[:,4]).min()))
 else:fit={k:np.nan for k in ['fit_t_start_ps','fit_t_end_ps','fit_samples','intensity_log_slope_ps_inv','field_sigma_fit_ps_inv','fit_R2','fit_inversion_min']}
 pred=np.log(np.sqrt(.81*peak)/seed)/sig
 fixed=cross(I,dt,1e-4);predfixed=np.log(.01/seed)/sig
 row=dict(case=name,source=str(folder.relative_to(ROOT)),energy_eV=e,detuning_eV=e-eth,step_ps=dt,terminal_ps=t[-1],maximum_I=peak,t_on_own81_ps=ton,fixed_I_1e_4_crossing_ps=fixed,fixed_I_1e_4_linear_prediction_ps=predfixed,t_on_own81_linear_prediction_ps=pred,nonlinear_criterion_correction_ps=ton-pred,sigma_ps_inv=sig,effective_seed_over_Esat=seed,instability_N=l['instability_N'],N_at_own81=float(np.interp(ton,tr[:,0],tr[:,4])),**fit)
 rows.append(row)
 for H in [15,20,30,40,60,80,200]:
  if H<=t[-1]+dt:
   J=I[:int(H/dt)+1];horizons.append(dict(case=name,energy_eV=e,horizon_ps=H,own_maximum_I=float(J.max()),own81_crossing_ps=cross(J,dt,.81*J.max())))
 if name.startswith('central_'):central.append(row)
 inputs.append(dict(case=name,parameters=par,metrics=met,source_files_sha256={n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ['parameters.tsv','metrics.tsv','trajectory.tsv']}))
 del I,t
write('trajectory_comparison.csv',rows);write('horizon_checks.csv',horizons);write('central_refinement.csv',sorted(central,key=lambda r:r['energy_eV']));(OUT/'inputs/selected_parameters.json').write_text(json.dumps(inputs,indent=2))
refinements=[]
for x in rows:
 match=None
 if x['case'].startswith('half_'):match='central_'+x['case'][5:]
 if x['case'].endswith('dt0.5'):match=x['case'].replace('dt0.5','dt1')
 if match:
  a=next(r for r in rows if r['case']==match);refinements.append(dict(case=x['case'],reference=match,delta_ton_ps=x['t_on_own81_ps']-a['t_on_own81_ps'],peak_relative_difference=x['maximum_I']/a['maximum_I']-1))
write('step_checks.csv',refinements)
# Separate observation-grid refinement from actual integration-step refinement.
cad=[]
for name in ['central_+0e+00','central_+1e-03','pair_left_dt1','pair_right_dt1']:
 tr,_,_=series[name];row=next(r for r in rows if r['case']==name)
 for stride in [1,5,20]:
  d=tr[::stride];lev=.81*row['maximum_I'];ii=np.flatnonzero((d[:-1,1]<lev)&(d[1:,1]>=lev));i=ii[0];crossing=d[i,0]+(lev-d[i,1])/(d[i+1,1]-d[i,1])*(d[i+1,0]-d[i,0]);cad.append(dict(case=name,output_stride=stride,effective_output_cadence_ps=float(np.median(np.diff(d[:,0]))),crossing_error_vs_native_ps=crossing-row['t_on_own81_ps']))
write('cadence_checks.csv',cad)
# Explicit carrier-coordinate check, no interpretation of imaginary shifts as growth.
ref=models[eth];cov=[];hbar=9.6/ref['k']
for e,m in models.items():
 A=m['A0']+m['D'];Ar=ref['A0']+ref['D'];delta=(e-eth)/hbar;res=A-Ar-1j*delta*np.eye(3)
 cov.append(dict(energy_eV=e,pure_carrier_shift_ps_inv=delta,noncovariant_matrix_residual_norm_ps_inv=float(np.linalg.norm(res)),sigma_ps_inv=lin[e]['sigma_ps_inv'],lab_mode_energy_eV=e-hbar*lin[e]['imag_ps_inv']))
write('carrier_check.csv',cov)
plt.rcParams.update({'font.size':10,'figure.constrained_layout.use':True})
def save(fig,name):
 fig.savefig(IMG/f'[OMEx]_{name}.pdf');fig.savefig(IMG/f'[OMEx]_{name}.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(10,7));es=sorted(lin)
ax[0,0].plot(es,[lin[e]['sigma_ps_inv'] for e in es],'.-');ax[0,0].set(ylabel='Amplitude growth σ (ps⁻¹)',xlabel='Probe energy (eV)')
ax[0,1].semilogy(es,[lin[e]['effective_output_seed_over_Esat'] for e in es],'.-');ax[0,1].set(ylabel='Effective dipolar seed / Esat',xlabel='Probe energy (eV)')
tracked=np.genfromtxt(OUT/'tables/tracked_modes.csv',names=True,delimiter=',')
for e,label in [(eth,'Exact threshold'),(2.5309774763480535,'Pair left'),(3.1201269537068566,'Pair right')]:
 a=tracked[tracked['energy_eV']==e];mode=int(lin[e]['mode']);a=a[a['mode']==mode];ax[1,0].plot(a['N'],a['sigma_ps_inv'],label=label)
ax[1,0].axhline(0,color='gray',lw=.7);ax[1,0].set(xlabel='Inversion N',ylabel='Tracked mode σ (ps⁻¹)');ax[1,0].legend()
rr=sorted([r for r in rows if not r['case'].startswith('half_') and not r['case'].endswith('dt0.5')],key=lambda r:r['energy_eV'])
ax[1,1].plot([r['energy_eV'] for r in rr],[r['t_on_own81_ps'] for r in rr],'o',label='Native own-81% crossing');ax[1,1].plot([r['energy_eV'] for r in rr],[r['t_on_own81_linear_prediction_ps'] for r in rr],'.-',label='Growing-mode extrapolation');ax[1,1].set(xlabel='Probe energy (eV)',ylabel='Time (ps)');ax[1,1].legend(fontsize=8);save(fig,'linear_and_delay')
fig,ax=plt.subplots(2,2,figsize=(10,7))
for name,label,color in [('pair_left_dt1','Left','#1565C0'),('pair_right_dt1','Right','#D97706')]:
 tr,J,dt=series[name];row=next(r for r in rows if r['case']==name);e=row['energy_eV'];s=lin[e];tx=np.arange(len(J))*dt
 ax[0,0].semilogy(tx,J,color=color,label=label);tf=np.linspace(1,28,500);ax[0,0].semilogy(tf,s['effective_output_seed_over_Esat']**2*np.exp(2*s['sigma_ps_inv']*tf),'--',color=color)
 ax[0,1].plot(tr[:,0],tr[:,4],color=color,label=label);ax[0,1].axhline(s['instability_N'],color=color,ls=':')
 ax[1,0].plot(tr[:,0],tr[:,1],color=color,label=label);ax[1,0].axvline(row['t_on_own81_ps'],color=color,ls=':')
 z=np.load(OUT/'tables'/f'linear_{e:.17g}.npz');a=tr[tr[:,0]<=2];ax[1,1].semilogy(a[:,0],a[:,1],color=color,label=label);ax[1,1].semilogy(z['t'],z['I'],'--',color=color,lw=.7)
ax[0,0].set(xlim=(0,28),ylim=(1e-19,.2),xlabel='Time (ps)',ylabel='I / Isat; dashed linear mode');ax[0,0].legend();ax[0,1].set(xlim=(0,28),xlabel='Time (ps)',ylabel='Inversion N');ax[1,0].set(xlim=(15,28),xlabel='Time (ps)',ylabel='I / Isat');ax[1,1].set(xlim=(0,2),ylim=(1e-20,1e-15),xlabel='Time (ps)',ylabel='Early I / Isat; dashed forced linear');save(fig,'pair_mechanisms')
fig,ax=plt.subplots(1,3,figsize=(12,3.7));cr=sorted(central,key=lambda x:x['energy_eV']);x=np.array([r['detuning_eV'] for r in cr]);y=np.array([r['t_on_own81_ps'] for r in cr]);center=next(r for r in cr if r['case']=='central_+0e+00')['t_on_own81_ps']
ax[0].plot(x*1e3,y,'o-');ax[0].set(xlabel='Detuning (meV)',ylabel='Native own-81% time (ps)')
mask=abs(x)<=1.1e-4;ax[1].plot(x[mask]*1e6,(y[mask]-center)*1e3,'o-');ax[1].set(xlabel='Detuning (µeV)',ylabel='Delay relative to exact threshold (fs)')
pairs=[]
for d in [.01,.001,.0001,.00001,.000001]:
 a=min(cr,key=lambda r:abs(r['detuning_eV']-d));b=min(cr,key=lambda r:abs(r['detuning_eV']+d));pairs.append(dict(abs_detuning_eV=d,even_delay_ps=(a['t_on_own81_ps']+b['t_on_own81_ps'])/2-center,odd_delay_ps=(a['t_on_own81_ps']-b['t_on_own81_ps'])/2))
write('central_even_odd.csv',pairs)
ax[2].loglog([r['abs_detuning_eV'] for r in pairs],[r['even_delay_ps'] for r in pairs],'o-',label='Even delay');ds=np.logspace(-6,-4,40);nuprime=(lin[eth+1e-6]['imag_ps_inv']-lin[eth-1e-6]['imag_ps_inv'])/2e-6;quad=nuprime**2/(2*lin[eth]['sigma_ps_inv']**3);ax[2].loglog(ds,quad*ds**2,'--',label='Linear-pole quadratic prediction');ax[2].set(xlabel='Absolute detuning (eV)',ylabel='Even delay increase (ps)');ax[2].legend();save(fig,'central_refinement')
fig,ax=plt.subplots(1,2,figsize=(10,3.7));edge=np.genfromtxt(OUT/'tables/edge_rates.csv',names=True,delimiter=',')
for a,(lo,hi) in zip(ax,[(2.451,2.456),(3.212,3.219)]):
 z=edge[(edge['energy_eV']>=lo)&(edge['energy_eV']<=hi)];a.plot(z['energy_eV'],z['sigma_ps_inv']);a.axhline(0,color='gray');a.set(xlabel='Probe energy (eV)',ylabel='Amplitude growth σ (ps⁻¹)')
save(fig,'edge_rates')
print('Analyzed',len(rows),'trajectories; saved tables and four figure sets')
