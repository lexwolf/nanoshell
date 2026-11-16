from pathlib import Path
import csv,json,numpy as np,re
paths=__import__("[OMEx]_paths");p=paths.DATA
def load(name):return list(csv.DictReader((p/'tables'/name).open()))
new=load('gain_dependence.csv');old=list(csv.DictReader((p/'sources/old_summaries/gain_dependence_linear_extended_summary.csv').open()));old={float(r['gain_ratio']):r for r in old}
rows=[]
for r in new:
 g=float(r['gain_ratio']);o=old[g];row={'gain_ratio':g}
 for n,k in [('peak_I','peak_Iem_over_Isat'),('peak_time_ps','peak_time_ps'),('fwhm_ps','fwhm_ps')]:
  row['old_'+n]=float(o[k]);row['new_'+n]=float(r[n]);row['difference_'+n]=float(r[n])-float(o[k])
 rows.append(row)
with (p/'tables/pulse_old_new.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('new width minimum',min(new,key=lambda r:float(r['fwhm_ps'])))
print('selected oldnew',[r for r in rows if r['gain_ratio'] in [1.05,1.09,1.1,1.15,1.2,1.25,1.5]])
probe=load('probe_scan.csv');print('probe weak',[next(r for r in probe if float(r['gain_ratio'])==g and int(r['index'])==0) for g in [1.05,1.1,1.15]])
freq=load('frequency_scan.csv');print('freqstates',[(g,sum(r['status']=='completed' for r in freq if float(r['gain_ratio'])==g),len([r for r in freq if float(r['gain_ratio'])==g])) for g in [1.05,1.1,1.15]])
# Preserve line-specific old text for editorial review, including all captions and quantitative prose.
files=[p/'sources/LaTeX/article/Gain-switched_emission_in_metal_nanoshell_spasers_OME.tex',*list((p/'sources/LaTeX/supplement').glob('*.tex')),p/'sources/LaTeX/Answer_to_Referee_Comments.tex']
lines=['# Source statements to review after regeneration','', 'Exact source excerpts below are an inventory, not proposed replacement text. Some remain valid; the comparison report records which numerical conclusions changed. Line numbers refer to the frozen sources.','']
for q in files:
 lines+=['## '+str(q.relative_to(p/'sources/LaTeX')),''];inside=False
 for i,line in enumerate(q.read_text().splitlines(),1):
  if '\\caption{' in line:inside=True
  relevant=(inside or (i>(100 if q.name=='Answer_to_Referee_Comments.tex' else 175) and re.search(r'(?:Fig\.|Table|10\^|[0-9]+\\,G|[0-9]+\\%|[0-9]+\\times|[0-9]+~ps|[0-9]+\\,\\mathrm\{ps\}|507|672|1\.09|2\.24|2\.58|12\.25|0\.90|1\.11|0\.317931|90\\%|81\\%)',line)))
  if relevant:lines += [f'- Line {i}: `{line}`']
  if inside and ('\\label' in line or '\\end{figure}' in line):inside=False
(p/'SOURCE_STATEMENT_INVENTORY.md').write_text('\n'.join(lines)+'\n')
