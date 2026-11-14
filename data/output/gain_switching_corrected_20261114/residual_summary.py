from pathlib import Path
import json,csv,collections,numpy as np
p=Path(__file__).resolve().parent;G=json.loads((p/'inputs/grids.json').read_text());rows=[];summary=[]
for g,grid in G.items():
 pairs=[]
 for i,e in enumerate(grid['frequencies']):
  rr=[]
  for tol in [1,.5]:
   name=f'residual_tight_g{g}_{i:04d}' if tol==.5 else (f'residual_default_g{g}_{i:04d}' if (p/f'runs/residual_default_g{g}_{i:04d}/result.json').exists() else f'freq_g{g}_{i:04d}')
   f=p/'runs'/name;r=json.loads((f/'result.json').read_text());meta=dict(l.split('\t',1) for l in (f/'diagnostic_metadata.tsv').read_text().splitlines());w=list(csv.DictReader((f/'diagnostic_windows.tsv').open(),delimiter='\t'));last=w[-1] if w else None
   v=float(last['relative_interval_range']) if last and last['relative_interval_range']!='null' else None
   row=dict(gain_ratio=float(g),frequency_eV=e,tolerance_factor=tol,run=name,run_status=r['metrics']['status'],diagnostic_status=meta['diagnostic_status'],diagnostic_complete=meta['diagnostic_complete']=='1',relative_interval_range=v,percent=100*v if v is not None else None,window_count=len(w),window_start_ps=last['start_ps'] if last else None,window_end_ps=last['end_ps'] if last else None,paired_stable=False)
   rr.append(row)
  good=all(x['diagnostic_complete'] and x['diagnostic_status']=='stable_interval_variation' and x['relative_interval_range'] is not None for x in rr)
  if good:
   va,vb=[x['relative_interval_range'] for x in rr];good=abs(va-vb)<=1e-8+.05*(va+vb)/2
  for x in rr:x['paired_stable']=good
  rows+=rr
  if good:pairs.append(max(rr,key=lambda r:r['relative_interval_range']))
 maximum=max(pairs,key=lambda r:r['relative_interval_range']) if pairs else None
 summary.append(dict(gain_ratio=float(g),recorded_energies=len(grid['frequencies']),paired_stable=len(pairs),unresolved_or_unstable=len(grid['frequencies'])-len(pairs),largest_paired_stable_percent=maximum['percent'] if maximum else None,energy_of_largest_eV=maximum['frequency_eV'] if maximum else None))
for name,rr in [('residual_selected_new.csv',rows),('supplement_table_S2.csv',summary)]:
 with (p/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
(p/'tables/residual_summary.json').write_text(json.dumps(summary,indent=2));print(summary)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(2,3,figsize=(12,5),sharex='col',gridspec_kw={'height_ratios':[5,1]})
for j,g in enumerate([1.05,1.1,1.15]):
 for tol,c in [(1,'#1565C0'),(.5,'#D97706')]:
  rr=[r for r in rows if r['gain_ratio']==g and r['tolerance_factor']==tol and r['paired_stable'] and r['percent']>0];axs[0,j].semilogy([r['frequency_eV'] for r in rr],[r['percent'] for r in rr],'.',ms=3,color=c,label=f'Tolerance ×{tol:g}')
 rr=[r for r in rows if r['gain_ratio']==g and r['tolerance_factor']==1];axs[1,j].scatter([r['frequency_eV'] for r in rr],[int(r['paired_stable']) for r in rr],s=3);axs[1,j].set_yticks([0,1],['Unresolved','Paired']);axs[1,j].set_xlabel('Excitation energy (eV)');axs[0,j].set_title(f'$G={g:.2f}G_{{th}}$');axs[0,j].legend(fontsize=8)
axs[0,0].set_ylabel('Late-window variation (%)');fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(p/'staging'/('residual_modulation_fig3.'+ext),dpi=180)
