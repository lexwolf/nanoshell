"""Index the reorganized review files without changing archived run provenance."""
from pathlib import Path
import hashlib,json,zipfile
paths=__import__('[OMEx]_paths')
p=paths.DATA
review=list(p.glob('*.md'))+[p/'RELOCATION.csv']
review += [q for q in (paths.ROOT/'scripts').rglob('*') if q.is_file() and '__pycache__' not in q.parts and (q.name.startswith(paths.PREFIX) or '[OMEx]_reference' in q.parts)]
review += [paths.DRIVER_SOURCE,paths.ROOT/'src/residual_diagnostic.hpp']
review += [q for q in paths.IMAGES.rglob('*') if q.is_file()]
for sub in ['inputs','sources/Data','sources/old_summaries','sources/nano_geo_matrix']:
 review += [q for q in (p/sub).rglob('*') if q.is_file()]
review += list((p/'sources/LaTeX').rglob('*.tex'))
review += list((p/'tables').glob('*.csv'))
review += [p/'tables'/n for n in ['coverage.json','validation.json','residual_summary.json','pdf_validation.json','plot_generation.json','gain_dependence_fits.json']]
review += [p/'compiler.txt',p/'build_command.json',p/'pilot/results.json',p/'pilot/PILOT.md',p/'REORGANIZATION_VALIDATION.json',p/'build_command_current.json',p/'compiler_current.txt']
review=sorted(set(q for q in review if q.is_file()))
rows=[]
for q in review:rows.append(hashlib.sha256(q.read_bytes()).hexdigest()+'  '+str(q.relative_to(paths.ROOT)))
(p/'REORGANIZED_SHA256SUMS').write_text('\n'.join(rows)+'\n')
with zipfile.ZipFile(p/(paths.PREFIX+'review_bundle.zip'),'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for q in review:z.write(q,q.relative_to(paths.ROOT))
 z.write(p/'REORGANIZED_SHA256SUMS',(p/'REORGANIZED_SHA256SUMS').relative_to(paths.ROOT))
print('Indexed',len(review),'review files; original SHA256SUMS remains the pre-reorganization record.')
