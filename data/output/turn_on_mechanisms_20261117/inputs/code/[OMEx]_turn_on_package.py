"""Snapshot investigation code and package compact review artifacts; keep native runs local."""
from pathlib import Path
import shutil,json,hashlib,subprocess,zipfile
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/turn_on_mechanisms_20261117';code=OUT/'inputs/code';code.mkdir(exist_ok=True)
files=list((ROOT/'scripts').glob('[[]OMEx]_turn_on_*'))+[ROOT/'src/[OMEx]_turn_on_matrix.cxx',ROOT/'src/[OMEx]_gain_switching.cxx',ROOT/'src/residual_diagnostic.hpp']
for f in files:
 if f.is_file():shutil.copy2(f,code/f.name)
shutil.copy2(ROOT/'docs/turn_on_mechanisms.md',OUT/'REPORT.md')
prov=json.loads((OUT/'inputs/provenance.json').read_text());prov['final_worktree_status']=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True);prov['analysis_sources_sha256']={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.is_file()};(OUT/'inputs/final_provenance.json').write_text(json.dumps(prov,indent=2))
imgs=ROOT/'img/output/[OMEx]_turn_on_mechanisms';keep=[f for f in OUT.rglob('*') if f.is_file() and f.name not in ['SHA256SUMS','[OMEx]_turn_on_review.zip']]+list(imgs.glob('*'))
with (OUT/'SHA256SUMS').open('w') as f:
 for q in sorted(keep):f.write(f'{hashlib.sha256(q.read_bytes()).hexdigest()}  {q.relative_to(ROOT)}\n')
review=[f for f in keep if f.suffix!='.bin' and ('/runs/' not in str(f) or f.name in ['parameters.tsv','metrics.tsv','command.json','result.json'])]
with zipfile.ZipFile(OUT/'[OMEx]_turn_on_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in review:z.write(f,str(f.relative_to(ROOT)))
 z.write(OUT/'SHA256SUMS',str((OUT/'SHA256SUMS').relative_to(ROOT)))
print('Packaged',len(review),'review artifacts; native histories retained separately.')
