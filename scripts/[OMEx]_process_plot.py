from pathlib import Path
import numpy as np,json,csv,gzip,collections,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
paths=__import__("[OMEx]_paths");p=paths.DATA;stage=paths.IMAGES;tables=p/'tables';diag=paths.DIAGNOSTICS
stage.mkdir(parents=True,exist_ok=True);diag.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':12,'axes.labelsize':13,'axes.linewidth':.9,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'ps.fonttype':42})
colors=['#0072B2','#D55E00','#009E73']; gains=[1.05,1.1,1.15];G=json.loads((p/'inputs/grids.json').read_text());pair=json.loads((p/'inputs/pair.json').read_text())
def read(q):return dict(line.split('\t',1) for line in q.read_text().splitlines())
def record(n):return json.loads((p/'runs'/n/'result.json').read_text())
def data(n):return np.loadtxt(p/'runs'/n/'trajectory.tsv')
def save(fig,name,folder=stage):
 fig.savefig(folder/(paths.PREFIX+name+'.pdf'),bbox_inches='tight');fig.savefig(folder/(paths.PREFIX+name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
def write(name,rows):
 with (tables/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));w.writeheader();w.writerows(rows)
def arrays(n):
 load=__import__("[OMEx]_native_io").load
 return load(p/'runs'/n)
# All original 500 normalized values, with newly computed raw E0 and all-step metrics.
probes=[]
for g in gains:
 for i,x in enumerate(G[str(g)]['probes']):
  r=record(f'probe_g{g}_{i:04d}');m=r['metrics'];a=r['parameters'];probes.append(dict(gain_ratio=g,index=i,E0_over_Esat=x,raw_E0=float(a['raw_E0']),Esat=float(a['Esat']),Isat=float(a['Isat']),t_on_ps=float(m['ton_81_peak_ps']),maximum_I=float(m['maximum_I']),status=m['status'],terminal_ps=float(m['terminal_ps'])))
write('probe_scan.csv',probes)
fig,ax=plt.subplots(figsize=(7,4.5))
for g,c in zip(gains,colors):
 rows=[r for r in probes if r['gain_ratio']==g];ax.semilogx([r['E0_over_Esat'] for r in rows],[r['t_on_ps'] for r in rows],color=c,label=f'{g:.2f} $G_{{th}}$')
ax.set(xlabel='$E_0/E_{sat}$',ylabel='$t_{on}$ (ps)',xlim=(10**-15.75,10**.75));ax.legend();fig.tight_layout();save(fig,'probe_middle_single')
fig,ax=plt.subplots(figsize=(7,4.5))
for g,c in zip(gains,colors):
 rows=[r for r in probes if r['gain_ratio']==g];ax.loglog([r['E0_over_Esat'] for r in rows],[r['maximum_I'] for r in rows],color=c,label=f'{g:.2f} $G_{{th}}$')
ax.set(xlabel='$E_0/E_{sat}$',ylabel='$I_{max}/I_{sat}$');ax.legend();fig.tight_layout();save(fig,'sm_probe_peak_response')
# Selected new residual default reruns override only same-parameter baseline entries.
freq=[]
for g in gains:
 mp=np.genfromtxt(p/f'runs/frequency_map_{g}/map.tsv',names=True,delimiter='\t')
 for i,e in enumerate(G[str(g)]['frequencies']):
  n=f'freq_g{g}_{i:04d}';other=f'residual_default_g{g}_{i:04d}'
  if (p/'runs'/other/'result.json').exists():n=other
  r=record(n);m=r['metrics'];valid=m.get('status')=='completed';freq.append(dict(gain_ratio=g,index=i,frequency_eV=e,emitted_analytic_eV=float(mp['emitted_eV'][i]),emitted_phase_eV=float(m.get('phase_emitted_eV',0)),t_on_ps=float(m['ton_81_peak_ps']) if valid else None,observed_crossing_ps=float(m['ton_81_peak_ps']),maximum_I=float(m['maximum_I']),status=m['status'],run=n,terminal_ps=float(m['terminal_ps'])))
write('frequency_scan.csv',freq)
for emitted,name in [(False,'freq_sweep_panels'),(True,'freq_sweep_LR_panels')]:
 fig,axs=plt.subplots(1,3,figsize=(11,4),sharey=True)
 for ax,g,c,letter in zip(axs,gains,colors,'abc'):
  rows=[r for r in freq if r['gain_ratio']==g]
  if not emitted:
   ax.plot([r['frequency_eV'] for r in rows],[r['t_on_ps'] if r['t_on_ps'] is not None else np.nan for r in rows],color=c)
  else:
   for left,ls in [(True,'-'),(False,'--')]:
    rr=[r for r in rows if (r['frequency_eV']<pair['threshold_eV'])==left];rr.sort(key=lambda r:r['emitted_phase_eV'])
    ax.plot([r['emitted_phase_eV'] for r in rr],[r['t_on_ps'] if r['t_on_ps'] is not None else np.nan for r in rr],ls,color=c,label='Left' if left else 'Right')
   ax.legend(fontsize=9);ax.ticklabel_format(useOffset=False,axis='x')
  missing=sum(r['t_on_ps'] is None for r in rows);ax.set_title(f'$G={g:.2f}G_{{th}}$');ax.set_ylim(0,100);ax.set_xlabel(r'$\hbar\omega_{em}$ (eV)' if emitted else r'$\hbar\omega$ (eV)');ax.text(.95,.08,f'({letter})',ha='right',transform=ax.transAxes);ax.text(.03,.94,f'{missing} incomplete / {len(rows)}',fontsize=8,transform=ax.transAxes,va='top');ax.locator_params(axis='x',nbins=4)
 axs[0].set_ylabel('$t_{on}$ (ps)');fig.tight_layout();save(fig,name)

def spectrum(x):
 x=np.asarray(x);r=np.where(x<2.1,1,np.where(x<2.4,-(x-2.4)/.3,np.where(x<2.8,0,np.where(x<3.3,(x-2.8)/.5,1))));g=np.where(x<1.9,0,np.where(x<2.1,(x-1.9)/.2,np.where(x<2.5,1,np.where(x<2.8,-(x-2.8)/.3,0))));b=np.where(x>2.5,1,np.where(x>2.4,-(x-2.4)/(-.1),0));f=np.where(x<1.5,0,np.where(x<1.8,.3+.7*(x-1.6)/.2,np.where(x<3.,1,np.where(x<3.41,.3+.7*(3.3-x)/.3,0))));return np.clip(np.stack([f*r,g,f*b],axis=-1),0,1)
spectral=ListedColormap(spectrum(np.linspace(1.8,3.6,1000)))
pulling_color=spectrum(pair['threshold_eV'])

# Autonomous stationary construction from a fresh library evaluation. No historical ordinates are used.
a=np.genfromtxt(p/'runs/pulling_map/map.tsv',names=True,delimiter='\t');scale=float(record('pulling_check_threshold')['parameters']['Isat']);I=np.maximum(0,a['I_raw'])/scale
fig=plt.figure(figsize=(7,4.5));gs=fig.add_gridspec(2,1,height_ratios=[8,1],hspace=.03);ax=fig.add_subplot(gs[0]);band=fig.add_subplot(gs[1],sharex=ax)
ax.plot(a['energy_eV'],I,color=pulling_color,lw=1.6);valid=(I>0)&(abs(a['Omega_normalized'])>0);em=a['emitted_eV'][valid];yy=I[valid];ix=np.argsort(em);ax.fill_between(em[ix],0,yy[ix],color=pulling_color,alpha=.4);ax.set_ylabel('$I_{em}/I_{sat}$');ax.set_xlim(1.8,3.6);ax.set_ylim(0,I.max()*1.12);ax.tick_params(labelbottom=False)
for level in I.max()*np.array([.15,.5,.78]):
 ids=np.flatnonzero((I[:-1]-level)*(I[1:]-level)<0)
 if len(ids)>=2:
  xs=[];targets=[]
  for j in [ids[0],ids[-1]]:
   w=(level-I[j])/(I[j+1]-I[j]);xs.append(a['energy_eV'][j]+w*(a['energy_eV'][j+1]-a['energy_eV'][j]));targets.append(a['emitted_eV'][j]+w*(a['emitted_eV'][j+1]-a['emitted_eV'][j]))
  for x,target in zip(xs,targets):ax.plot([x,x],[0,level],color=spectrum(x),lw=1);ax.annotate('',xy=(target,level),xytext=(x,level),arrowprops=dict(arrowstyle='->',color=spectrum(x)))
band.imshow(np.arange(1000)[None,:],extent=[1.8,3.6,0,1],aspect='auto',cmap=spectral);band.set_yticks([]);band.set_xlabel('Energy (eV)');save(fig,'pulling')
np.savetxt(tables/'pulling_corrected.tsv',np.c_[a['energy_eV'],a['emitted_eV'],I,a['Nss']],header='excitation_eV emitted_eV I_over_Isat Nss')
# Paired actual excitation energies, using analytic pairing and dynamical checks separately.
fig,ax=plt.subplots(figsize=(7,4.5));pairrows=[]
for side,c,ls in [('left','#1565C0','-'),('right','#D97706','--')]:
 n=f'pair_{side}_dt1';r=record(n);d=data(n);m=r['metrics'];ton=float(m['ton_81_peak_ps']);peak=float(m['maximum_I']);mask=d[:,0]<=40
 ax.plot(d[mask,0],d[mask,1]*1000,ls,color=c,lw=1.8,label=f'{side.capitalize()}: {pair[side+"_eV"]:.10f} eV  ($t_{{on}}$={ton:.2f} ps)');ax.axhline(.81*peak*1000,color='#888888',lw=.6,ls=':');ax.plot([ton,ton],[0,.81*peak*1000],':',color=c);ax.plot(ton,.81*peak*1000,'o',color=c,ms=4)
 final=float(d[-int(20/.000342818):,1].mean());pairrows.append(dict(branch=side,energy_eV=pair[side+'_eV'],energy_ratio=pair[side+'_eV']/pair['threshold_eV'],emitted_phase_eV=float(m['phase_emitted_eV']),final_20ps_I=final,final_last100000_saved_I=float(d[-100000:,1].mean()),t_on_ps=ton,peak=peak))
ax.set(xlim=(0,40),ylim=(0,1.15*max(r['peak'] for r in pairrows)*1000),xlabel='$t$ (ps)',ylabel='$10^3 I_{em}/I_{sat}$',title='$G=1.10G_{th}$');ax.legend(fontsize=9,loc='lower right');fig.tight_layout();save(fig,'paired_frequency_temporal_traces');write('paired_comparison.csv',pairrows)
# Pulse metrics from every integration step, not display samples.
pulserows=[];cyclerows=[]
for q in sorted((p/'runs').glob('pulse_g*')):
 g=float(q.name.split('g')[-1]);r=record(q.name);cyc=np.genfromtxt(q/'cycles.tsv',names=True,delimiter='\t');last=cyc[-1]
 pulserows.append(dict(gain_ratio=g,peak_I=float(last['peak_I']),peak_time_ps=float(last['peak_time_ps']),fwhm_ps=float(last['fwhm_ps']),half_rise_ps=float(last['half_rise_ps']),half_fall_ps=float(last['half_fall_ps']),Esat=float(r['parameters']['Esat']),raw_E0=float(r['parameters']['raw_E0']),status=r['metrics']['status']))
 for row in cyc:cyclerows.append(dict(gain_ratio=g,**{k:float(row[k]) for k in cyc.dtype.names}))
write('gain_dependence.csv',pulserows);write('all_cycle_metrics.csv',cyclerows)
fig,axs=plt.subplots(2,1,figsize=(7,6),sharex=True);d=data('pulse_g1.10');m=d[:,0]<=100
axs[0].plot(d[m,0],d[m,4],color='#8436A8',label='$N$');axs[0].plot(d[m,0],d[m,5],'--',color='#D97706',label=r'$\widetilde N$');axs[0].axhline(1/1.1,color='gray',ls=':');axs[0].legend();axs[0].set_ylabel('Inversion');axs[1].plot(d[m,0],d[m,1],color='#1565C0');axs[1].set(xlabel='$t$ (ps)',ylabel='$I_{em}/I_{sat}$');
# First upward and subsequent downward actual-N threshold crossings.
ix=np.flatnonzero((d[:-1,4]-1/1.1)*(d[1:,4]-1/1.1)<0)
for j in ix[:2]:
 t=d[j,0]+(1/1.1-d[j,4])*(d[j+1,0]-d[j,0])/(d[j+1,4]-d[j,4]);
 for ax in axs:ax.axvline(t,color='gray',ls='--',lw=.7)
for ax,l in zip(axs,'ab'):ax.text(.02,.88,f'({l})',transform=ax.transAxes)
fig.tight_layout();save(fig,'gain_switching_mechanism_1.1Gth')
fig,axs=plt.subplots(2,3,figsize=(11,6),sharex=True)
for j,g in enumerate(gains):
 d=data(f'pulse_g{g:.2f}');d=d[d[:,0]<50];axs[0,j].plot(d[:,0],d[:,4],color='#8436A8');axs[0,j].plot(d[:,0],d[:,5],'--',color='#D97706');axs[0,j].axhline(1/g,color='gray',ls=':');axs[0,j].set_title(f'$G_{{max}}={g:.2f}G_{{th}}$');axs[1,j].plot(d[:,0],d[:,1],color=colors[j]);axs[1,j].set_xlabel('$t$ (ps)');axs[1,j].ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
 for row in [0,1]:axs[row,j].text(.04,.88,f'({chr(97+row*3+j)})',transform=axs[row,j].transAxes)
axs[0,0].set_ylabel('Inversion');axs[1,0].set_ylabel('$I_{em}/I_{sat}$');fig.tight_layout();save(fig,'gain_switching_mechanism_composite')
fig,axs=plt.subplots(1,3,figsize=(11,5));xx=[r['gain_ratio'] for r in pulserows]
for ax,col,l,c,letter in zip(axs,['fwhm_ps','peak_I','peak_time_ps'],['FWHM (ps)','$I_{peak}/I_{sat}$','Peak time (ps)'],['#1565C0','#D97706','#8436A8'],'abc'):
 ax.plot(xx,[r[col] if r[col]>=0 else np.nan for r in pulserows],'-o',color=c,ms=3);ax.set(xlabel='$G_{max}/G_{th}$',ylabel=l,xlim=(1.04,1.51));ax.text(.04,.93,f'({letter})',transform=ax.transAxes)
axs[1].set_ylim(bottom=0);ins=axs[1].inset_axes([.10,.48,.55,.39]);ins.semilogy(xx,[r['peak_I'] for r in pulserows],'-o',color='#D97706',ms=2);ins.set_xlim(1.045,1.18);ins.set_ylim(min(r['peak_I'] for r in pulserows)*.5,max(r['peak_I'] for r in pulserows if r['gain_ratio']<=1.18)*2);ins.tick_params(labelsize=7);ins.set_title('Weak peaks (log scale)',fontsize=8);ins.grid(alpha=.2)
fig.suptitle(r'$\omega=\omega_{th},\ T=50$ ps, $E_0/E_{sat}=10^{-9}$');fig.tight_layout();save(fig,'gain_dependence_linear_extended')
# Supplement S2 and S3 retain their 3x2 arrangement and local-time windows.
for overlay,name in [(False,'sm_gain_switching_transition_mechanism'),(True,'sm_gain_switching_transition_cycles')]:
 fig,axs=plt.subplots(3,2,figsize=(9,9),sharex=True)
 for j,g in enumerate([1.15,1.2,1.25]):
  d=data(f'pulse_g{g:.2f}')
  for cycle in (range(4) if overlay else [3]):
   m=(d[:,0]>=cycle*50+12)&(d[:,0]<=cycle*50+38);a=d[m];x=a[:,0]-cycle*50;ls=['-','--','-.',':'][cycle] if overlay else '-';axs[j,0].plot(x,a[:,4],ls,lw=1,label=f'Cycle {cycle+1}');axs[j,1].plot(x,a[:,1],ls,lw=1,label=f'Cycle {cycle+1}')
  axs[j,0].plot(x,.5*(1-np.cos(2*np.pi*x/50)),'--',color='#D97706',label='Target');axs[j,0].axhline(1/g,color='gray',ls=':');axs[j,0].set_ylabel(f'$G_{{max}}/G_{{th}}={g:.2f}$\nInversion');axs[j,1].set_ylabel('$I_{em}/I_{sat}$');axs[j,1].axvline(25,color='gray',lw=.6);axs[j,1].ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
  if not overlay:
   r=next(r for r in pulserows if r['gain_ratio']==g);axs[j,1].plot(r['peak_time_ps'],r['peak_I'],'o',color='#1565C0',ms=4);axs[j,1].plot([r['half_rise_ps'],r['half_fall_ps']],[r['peak_I']/2]*2,color='#1565C0',lw=2)
 axs[0,0].legend(fontsize=8);axs[-1,0].set_xlabel('$s=t-kT$ (ps)');axs[-1,1].set_xlabel('$s=t-kT$ (ps)');axs[-1,0].set_xlim(12,38);fig.tight_layout();save(fig,name)
write('supplement_table_S1.csv',[r for r in pulserows if r['gain_ratio'] in [1.15,1.2,1.25]])
repeat=[]
for g in [1.15,1.2,1.25]:
 rr=[r for r in cyclerows if r['gain_ratio']==g];repeat.append(dict(gain_ratio=g,peak_time_span_ps=np.ptp([r['peak_time_ps'] for r in rr]),FWHM_span_ps=np.ptp([r['fwhm_ps'] for r in rr]),peak_relative_span_percent=100*np.ptp([r['peak_I'] for r in rr])/rr[-1]['peak_I']))
write('cycle_repeatability.csv',repeat)
# Reference refinement and cadence evidence; retain full native intensity arrays separately.
checks=[]
for side in ['left','right']:
 a=arrays(f'pair_{side}_dt1');b=arrays(f'pair_{side}_dt0.5');step=float(record(f'pair_{side}_dt1')['parameters']['integration_step_ps']);t=np.arange(len(a))*step;m=(t>=18)&(t<=22);diff=np.abs(a-b[::2][:len(a)])
 checks.append(dict(case=side,max_abs_I_difference_18_22=float(diff[m].max()),relative_to_peak=float(diff[m].max()/a.max()),ton_dt1=record(f'pair_{side}_dt1')['metrics']['ton_81_peak_ps'],ton_dt05=record(f'pair_{side}_dt0.5')['metrics']['ton_81_peak_ps'],phase_energy_dt1=record(f'pair_{side}_dt1')['metrics']['phase_emitted_eV'],phase_energy_dt05=record(f'pair_{side}_dt0.5')['metrics']['phase_emitted_eV']))
a=arrays('pair_left_dt1');b=arrays('pair_left_dense');(diag/'cadence_check.json').write_text(json.dumps({'native_I_bit_identical':bool(np.array_equal(a,b))},indent=2));write('pair_refinement.csv',checks)
for g in [1.05,1.1,1.15,1.2,1.25,1.5]:
 a=np.genfromtxt(p/f'runs/pulse_g{g:.2f}/cycles.tsv',names=True,delimiter='\t')[-1];b=np.genfromtxt(p/f'runs/pulse_half_g{g:.2f}/cycles.tsv',names=True,delimiter='\t')[-1];checks.append(dict(case=f'pulse_{g}',peak_relative_difference=float(abs(a['peak_I']-b['peak_I'])/a['peak_I']),FWHM_abs_difference_ps=float(abs(a['fwhm_ps']-b['fwhm_ps'])),peak_time_abs_difference_ps=float(abs(a['peak_time_ps']-b['peak_time_ps']))))
write('refinement_checks.csv',checks)
fig,axs=plt.subplots(4,1,figsize=(9,10),sharex=True)
for n in ['pair_left_dt1','pair_left_dt0.5','pair_left_dense']:
 a=data(n);a=a[(a[:,0]>=19.98)&(a[:,0]<=20.02)]
 for ax,y in zip(axs,[a[:,1],a[:,2],a[:,3],a[:,4]]):ax.plot(a[:,0],y,'.-',ms=2,lw=.8,label=n)
for ax,l in zip(axs,['$I/I_{sat}$','Re(dip)','Im(dip)','N']):ax.set_ylabel(l);ax.grid(alpha=.2)
axs[0].legend(fontsize=8);axs[-1].set_xlabel('$t$ (ps)');fig.tight_layout();save(fig,'fig5_resolved_feature',diag)
# Historical point-marker close-up is a comparison only, never a source for corrected figures.
hist=p/'sources/historical_pair/0.9omegath.dat'
if hist.exists():
 h=np.loadtxt(hist);fig,axs=plt.subplots(2,1,figsize=(9,6))
 for ax,lim in zip(axs,[(18,22),(20,20.01)]):
  m=(h[:,0]>=lim[0])&(h[:,0]<=lim[1]);a=h[m];ax.plot(a[:,0],(a[:,3]**2+a[:,4]**2)/.317931**2,'.-',ms=2,lw=.7,label='Historical: old seed and old divisor');ax.set_ylabel('Historical $I/I_{sat}$');ax.grid(alpha=.2);ax.legend(fontsize=8)
 axs[-1].set_xlabel('$t$ (ps)');fig.tight_layout();save(fig,'fig5_historical_samples',diag)
fig,ax=plt.subplots(figsize=(8,4))
for n in ['pair_left_dt1','pair_left_dt0.5','pair_left_dense']:
 a=data(n);a=a[(a[:,0]>=18)&(a[:,0]<=22)];ax.plot(a[:,0],a[:,1],'.-',ms=1,lw=.8,label=n)
ax.set(xlabel='$t$ (ps)',ylabel='$I/I_{sat}$');ax.legend(fontsize=8);fig.tight_layout();save(fig,'fig5_18_22ps',diag)
(p/'tables/plot_generation.json').write_text(json.dumps({'figures':sorted(f.name for f in stage.glob('*.pdf')),'metric_resolution':'every_integration_step_native_f64','plot_source':'new trajectories and new autonomous library construction','note':'Figs3/4 use identical new run selection; Fig4 abscissa is measured phase rotation of the same trajectory at every integration step over its final complete 10 ps block plus terminal partial block; analytic reference is exported separately; no historical emission coordinates reused.'},indent=2))
print('PLOTS_COMPLETE',len(list(stage.glob('*.pdf'))),flush=True)
