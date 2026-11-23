"""Native-step comparison of separate full-metal solver with saved first-order records."""
from pathlib import Path
import csv,json,gzip,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
m=__import__('[OMEx]_metal_envelope_validation');ROOT=m.ROOT;OUT=m.OUT;IMG=m.IMG;IMG.mkdir(exist_ok=True);c=m.constants();eth=c['old_threshold_eV'];ep=c['plasma_eV'];k=c['k_ps_inv'];G=1.1*c['old_Gth']
oldbase=ROOT/'data/output/gain_switching_corrected_20261114';prior=ROOT/'data/output/turn_on_mechanisms_20261117'
loadold=__import__('[OMEx]_native_io').load
read=lambda f:dict(x.split('\t',1) for x in f.read_text().splitlines())
def load(f):
 with gzip.open(f/'intensity.f64.gz','rb') as z:return np.frombuffer(z.read(),dtype='<f8')
def crossing(I,h,lev):
 ix=np.flatnonzero((I[:-1]<lev)&(I[1:]>=lev))
 if not len(ix):return np.nan
 i=ix[0];return (i+(lev-I[i])/(I[i+1]-I[i]))*h
rows=[];horizons=[];cache={}
for f in sorted((OUT/'runs').glob('*')):
 if not (f/'metrics.json').exists() or f.name.startswith('passive_'):continue
 a=json.loads((f/'metrics.json').read_text());tr=np.loadtxt(f/'trajectory.tsv');I=load(f);dt=a['step_ps'];peak=I.max();ton=crossing(I,dt,.81*peak);assert abs(ton-a['ton81_ps'])<1e-10
 tail=tr[tr[:,0]>=tr[-1,0]-5];dip=tail[:,2]+1j*tail[:,3];ph=np.unwrap(np.angle(dip));rate=np.polyfit(tail[:,0],ph,1)[0];em=a['reference_eV']-ep/k*rate
 z=I[int(np.ceil((a['terminal_ps']-5)/dt)):];late=z.mean();t=np.arange(len(I))*dt;mask=(I>=1e-6)&(I<=1e-5)&(t>=1);fit=np.polyfit(t[mask],np.log(I[mask]),1) if mask.sum()>20 else [np.nan,np.nan]
 lam=eigvals(m.matrix(1,a['G'],a['gain_center_eV']/ep,a['reference_eV']/ep)[0]);sigma=max(lam.real)*k
 row=dict(case=f.name,formulation='second_order',probe_eV=a['probe_eV'],gain_center_eV=a['gain_center_eV'],G=a['G'],dt_normalized=a['dt_normalized'],observation_ps=a['terminal_ps'],preparation=a['preparation'],sigma_ps_inv=sigma,fitted_amplitude_sigma_ps_inv=fit[0]/2,maximum_J=peak,ton81_ps=ton,fixed_J_1e_4_ps=a['fixed_J_1e_4_ps'],late_mean_J=late,late_relative_ptp=np.ptp(z)/late,late_mean_N=tail[:,4].mean(),emitted_eV=em,legacy_time_label=False)
 rows.append(row)
 stride=max(1,int(.001/dt));cache[f.name]=(tr,np.arange(0,len(I),stride)*dt,I[::stride].copy(),row)
 for H in [20,25,30]:
  J=I[:int(H/dt)+1];horizons.append(dict(case=f.name,H_ps=H,maximum_J=J.max(),ton81_ps=crossing(J,dt,.81*J.max())))
 del I,t
# Same 30 ps intensity reference for the direct comparison; late state also checked on original full history.
selected=[('threshold',prior/'runs/central_+0e+00'),('left',oldbase/'runs/pair_left_dt1'),('right',oldbase/'runs/pair_right_dt1')]
for d in [-.01,-.001,-.0001,-.00001,.00001,.0001,.001,.01]:selected.append((f'central_{d:+.0e}',prior/'runs'/f'central_{d:+.0e}'))
oldlin={float(x['energy_eV']):x for x in csv.DictReader((prior/'tables/linear_summary.csv').open())}
for tag,f in selected:
 p=read(f/'parameters.tsv');met=read(f/'metrics.tsv');e=float(p['frequency_eV']);dt=float(p['integration_step_ps']);allI=loadold(f);I=allI[:int(30/dt)+1];tr=np.loadtxt(f/'trajectory.tsv');tail=tr[tr[:,0]>=tr[-1,0]-5];d=tail[:,2]+1j*tail[:,3];ph=np.unwrap(np.angle(d));em=e-ep/k*np.polyfit(tail[:,0],ph,1)[0];z=allI[int((len(allI)*dt-5)/dt):]
 row=dict(case='old_'+tag,formulation='first_order',probe_eV=e,gain_center_eV=eth,G=float(p['Gmax']),dt_normalized=float(p['integration_step_normalized']),observation_ps=(len(I)-1)*dt,preparation='reduced_q_zero_slaved_derivative',sigma_ps_inv=float(oldlin[e]['sigma_ps_inv']),fitted_amplitude_sigma_ps_inv=np.nan,maximum_J=I.max(),ton81_ps=crossing(I,dt,.81*I.max()),fixed_J_1e_4_ps=crossing(I,dt,1e-4),late_mean_J=z.mean(),late_relative_ptp=np.ptp(z)/z.mean(),late_mean_N=tail[:,4].mean(),emitted_eV=em,legacy_time_label=True)
 rows.append(row);stride=max(1,int(.001/dt));cache['old_'+tag]=(tr,np.arange(0,len(I),stride)*dt,I[::stride].copy(),row)
 del allI,I
m.write('trajectory_comparison.csv',rows);m.write('horizon_checks.csv',horizons)
# Static thresholds and spectra at exactly the same physical parameters and probes.
oldsum={float(x['energy_eV']):x for x in csv.DictReader((prior/'tables/linear_summary.csv').open())};static=[];modes=[];fullth=json.loads((OUT/'tables/thresholds.json').read_text())
for e in sorted(set(x['probe_eV'] for x in rows if not x['case'].startswith('matched_'))):
 if e not in oldsum:continue
 vals=eigvals(m.matrix(1,G,eth/ep,e/ep)[0]);dom=vals[vals.real.argmax()]
 for N in [0,.5,fullth['critical_N'],1]:
  vs=eigvals(m.matrix(N,G,eth/ep,e/ep)[0])
  for j,v in enumerate(sorted(vs,key=lambda x:x.imag)):modes.append(dict(probe_reference_eV=e,N=N,mode_imag_order=j,sigma_ps_inv=v.real*k,nu_ps_inv=v.imag*k,lab_energy_eV=e-ep*v.imag))
 o=oldsum[e];static.append(dict(probe_reference_eV=e,old_sigma_ps_inv=float(o['sigma_ps_inv']),new_sigma_ps_inv=dom.real*k,old_critical_N=float(o['instability_N']),new_critical_N=fullth['critical_N'],old_required_G=G*float(o['instability_N']),new_required_G=fullth['fixed_gain_center_Gth']))
m.write('spectral_comparison.csv',static);m.write('all_full_modes.csv',modes)
# Quantify the central even delay at successively smaller detunings.
center=cache['threshold_dt0.125'][3]['ton81_ps'];centerold=cache['old_threshold'][3]['ton81_ps'];smooth=[]
for d in [.01,.001,.0001,.00001]:
 a=cache[f'central_{d:+.0e}'][3];b=cache[f'central_{-d:+.0e}'][3];oa=cache[f'old_central_{d:+.0e}'][3];ob=cache[f'old_central_{-d:+.0e}'][3]
 smooth.append(dict(detuning_eV=d,new_even_delay_ps=(a['ton81_ps']+b['ton81_ps'])/2-center,old_even_delay_ps=(oa['ton81_ps']+ob['ton81_ps'])/2-centerold,new_odd_delay_ps=(a['ton81_ps']-b['ton81_ps'])/2))
m.write('central_smoothness.csv',smooth)
comparisons=[]
for tag in ['threshold','left','right']:
 a=cache['old_'+tag][3];b=cache[tag+'_dt0.125'][3];s=cache['slaved_'+tag][3]
 comparisons.append(dict(case=tag,old_ton81_30ps=a['ton81_ps'],new_rest_ton81_30ps=b['ton81_ps'],new_slaved_ton81_30ps=s['ton81_ps'],new_minus_old_ps=b['ton81_ps']-a['ton81_ps'],preparation_control_shift_ps=s['ton81_ps']-b['ton81_ps'],old_late_J=a['late_mean_J'],new_late_J=b['late_mean_J'],old_emitted_eV=a['emitted_eV'],new_emitted_eV=b['emitted_eV']))
m.write('key_comparison.csv',comparisons)
pairsummary={'old_pair_delay_ps':cache['old_right'][3]['ton81_ps']-cache['old_left'][3]['ton81_ps'],'new_same_pair_rest_delay_ps':cache['right_dt0.125'][3]['ton81_ps']-cache['left_dt0.125'][3]['ton81_ps'],'new_same_pair_slaved_delay_ps':cache['slaved_right'][3]['ton81_ps']-cache['slaved_left'][3]['ton81_ps'],'new_matched_symmetric_pair_delay_ps':cache['matched_right'][3]['ton81_ps']-cache['matched_left'][3]['ton81_ps'],'autonomous_map':'Carrier-invariant constant emission frequency at the selected gain center and G. No unique folded left/right pairing remains.'}
(OUT/'tables/pair_summary.json').write_text(json.dumps(pairsummary,indent=2))
plt.rcParams.update({'font.size':10,'figure.constrained_layout.use':True})
def save(fig,name):
 for ext in ['pdf','png']:fig.savefig(IMG/f'[OMEx]_{name}.{ext}',dpi=180)
 plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(10,7))
for e in ['left','right']:
 for old,ls in [(True,'--'),(False,'-')]:
  name='old_'+e if old else e+'_dt0.125';tr,t,I,r=cache[name];col='#1565C0' if e=='left' else '#D97706';ax[0,0].plot(t,I,ls,color=col,label=('Reduced ' if old else 'Full ')+e);ax[0,1].plot(tr[:,0],tr[:,4],ls,color=col)
a=sorted(static,key=lambda r:r['probe_reference_eV']);xx=[r['probe_reference_eV'] for r in a];ax[1,0].plot(xx,[r['old_sigma_ps_inv'] for r in a],'o--',label='Reduced metal');ax[1,0].plot(xx,[r['new_sigma_ps_inv'] for r in a],'.-',label='Full metal')
for old,ls in [(True,'--'),(False,'-')]:
 rr=[]
 for d in [-.01,-.001,-.0001,-.00001,0,.00001,.0001,.001,.01]:
  name=('old_threshold' if old else 'threshold_dt0.125') if d==0 else ('old_' if old else '')+f'central_{d:+.0e}';rr.append((d*1000,cache[name][3]['ton81_ps']))
 ax[1,1].plot(*zip(*rr),'o'+ls,label='Reduced' if old else 'Full')
ax[0,0].set(xlim=(8,28),xlabel='Time (ps)',ylabel='J = I/Isat');ax[0,0].legend(fontsize=8);ax[0,1].set(xlim=(8,28),ylim=(.88,1.01),xlabel='Time (ps)',ylabel='Inversion N');ax[1,0].set(xlabel='Probe / reference energy (eV)',ylabel='Full-inversion amplitude growth (ps⁻¹)');ax[1,0].legend();ax[1,1].set(xlabel='Detuning from old threshold (meV)',ylabel='Own-81% turn-on (ps)');ax[1,1].legend();save(fig,'metal_comparison')
fig,ax=plt.subplots(1,3,figsize=(12,3.6))
for key,mark in [('old_even_delay_ps','s'),('new_even_delay_ps','o')]:ax[0].loglog([x['detuning_eV'] for x in smooth],[x[key] for x in smooth],mark+'-',label='Reduced' if key.startswith('old') else 'Full')
ax[0].set(xlabel='Absolute detuning (eV)',ylabel='Even delay increase (ps)');ax[0].legend()
for tag,col in [('left','#1565C0'),('right','#D97706')]:
 for prep,ls in [('rest','-'),('slaved','--')]:
  name=tag+'_dt0.125' if prep=='rest' else 'slaved_'+tag;tr,t,I,r=cache[name];ax[1].plot(t,I,ls,color=col,label=tag+' '+prep)
ax[1].set(xlim=(11,15),xlabel='Time (ps)',ylabel='J = I/Isat');ax[1].legend(fontsize=8)
for tag,col in [('left','#1565C0'),('right','#D97706')]:
 tr,t,I,r=cache['matched_'+tag];ax[2].plot(t,I,color=col,label='Matched '+tag)
ax[2].set(xlim=(12,14),xlabel='Time (ps)',ylabel='Matched-gain J = I/Isat');ax[2].legend(fontsize=8);save(fig,'metal_preparation_and_central')
print(json.dumps(comparisons,indent=2));print(json.dumps(pairsummary,indent=2))
