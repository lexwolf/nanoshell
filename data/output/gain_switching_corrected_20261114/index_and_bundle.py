from pathlib import Path
import hashlib,json,zipfile,subprocess
p=Path(__file__).resolve().parent
# Native histories have independently verified uncompressed SHA-256s per run.
# Hash the other scientific records, including output trajectories, once.
rows=[]
for q in sorted(p.rglob('*')):
 if not q.is_file() or q.name in ['SHA256SUMS','review_bundle.zip'] or q.name in ['intensity.sf64.gz','intensity.f64','intensity.f64.gz'] or '__pycache__' in q.parts:continue
 h=hashlib.sha256()
 with q.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 rows.append(h.hexdigest()+'  '+str(q.relative_to(p)))
(p/'SHA256SUMS').write_text('\n'.join(rows)+'\n')
review=[]
review.extend(p.glob('*.md'));review.extend(p.glob('*.py'));review.extend(p.glob('*.cxx'));review.append(p/'residual_diagnostic.hpp')
review.extend((p/'staging').glob('*'));review.extend((p/'diagnostics').glob('*.pdf'));review.extend((p/'diagnostics').glob('*.png'))
review.extend((p/'tables').glob('*.csv'));review.extend(q for q in (p/'sources/Data').rglob('*') if q.is_file())
review.extend(q for q in (p/'sources/old_summaries').glob('*') if q.is_file())
review.extend((p/'sources/LaTeX').rglob('*.tex'))
review.extend(p/'tables'/n for n in ['coverage.json','validation.json','residual_summary.json','pdf_validation.json','plot_generation.json'])
review.extend((p/'inputs').glob('*'));review.extend(q for q in (p/'sources/nano_geo_matrix').rglob('*') if q.is_file());review.append(p/'build_command.json');review.append(p/'compiler.txt');review.append(p/'pilot/results.json');review.append(p/'pilot/PILOT.md')
with zipfile.ZipFile(p/'review_bundle.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for q in sorted(set(review)):
  if q.is_file():z.write(q,q.relative_to(p))
print('Indexed records:',len(rows),'review bundle bytes:',(p/'review_bundle.zip').stat().st_size)
for repo in ['.', 'extern/nano_geo_matrix']:
 print(repo,subprocess.check_output(['git','status','--short'],cwd=p.parents[2]/repo,text=True).strip() or 'clean')
