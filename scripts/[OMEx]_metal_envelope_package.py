"""Snapshot sources/provenance and build a compact review bundle, excluding native histories."""
from pathlib import Path
import json,hashlib,subprocess,shutil,zipfile,ast
m=__import__('[OMEx]_metal_envelope_validation');ROOT=m.ROOT;OUT=m.OUT;IMG=m.IMG
assert json.loads((OUT/'tables/final_validation.json').read_text())['passed']
code=OUT/'inputs/code';code.mkdir(exist_ok=True)
files=list((ROOT/'scripts').glob('[[]OMEx]_metal_envelope_*.py'))+[ROOT/'src/[OMEx]_metal_envelope_comparison.cxx',ROOT/'src/[OMEx]_gain_switching.cxx',ROOT/'src/residual_diagnostic.hpp']
for q in files:
 shutil.copy2(q,code/q.name)
 if q.suffix=='.py':ast.parse(q.read_text())
shutil.copy2(ROOT/'docs/metal_envelope_validation.md',OUT/'REPORT.md')
# Reference calibration: authoritative shared-scale and ordinary-silver implementation.
for q in [ROOT/'extern/nano_geo_matrix/docs/saturation_scale.md',ROOT/'extern/nano_geo_matrix/include/nano_geo_matrix/core/saturation_scale.hpp']:
 shutil.copy2(q,OUT/'inputs'/q.name)
base=Path('/home/alessandro/Nextcloud/Work/Research/Papers/Active/Gain Switching/round 1/LaTeX')
for rel in ['Answer_to_Referee_Comments.tex','article/Gain-switched_emission_in_metal_nanoshell_spasers_OME.tex']:
 shutil.copy2(base/rel,OUT/'inputs'/Path(rel).name)
# Preserve exactly which first-order sources were read, without duplicating their large histories.
old=[ROOT/'data/output/gain_switching_corrected_20261114/runs'/n for n in ['pair_left_dt1','pair_right_dt1','pair_left_dt0.5','pair_right_dt0.5']]
old += [ROOT/'data/output/turn_on_mechanisms_20261117/runs'/n for n in ['central_+0e+00','half_+0e+00']+[f'central_{d:+.0e}' for d in [-.01,-.001,-.0001,-.00001,.00001,.0001,.001,.01]]]
reused=[]
for f in old:
 reused.append({'folder':str(f.relative_to(ROOT)),'files_sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in f.iterdir() if q.is_file() and (q.name in ['parameters.tsv','metrics.tsv','trajectory.tsv'] or q.name.startswith('intensity.'))}})
(OUT/'inputs/reused_sources.json').write_text(json.dumps(reused,indent=2))
prov={'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'submodule_HEAD':subprocess.check_output(['git','-C','extern/nano_geo_matrix','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'submodule_status':subprocess.check_output(['git','-C','extern/nano_geo_matrix','status','--short'],cwd=ROOT,text=True),'worktree_status':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),'compiler':subprocess.check_output(['g++','--version'],text=True),'source_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in files},'binary_sha256':hashlib.sha256((OUT/'comparison.bin').read_bytes()).hexdigest(),'new_runs':sum((f/'result.json').exists() for f in (OUT/'runs').glob('*'))}
assert prov['submodule_status']==''
(OUT/'inputs/provenance.json').write_text(json.dumps(prov,indent=2))
keep=[q for q in OUT.rglob('*') if q.is_file() and q.name not in ['SHA256SUMS','[OMEx]_metal_envelope_review.zip']]+list(IMG.glob('*'))
with (OUT/'SHA256SUMS').open('w') as f:
 for q in sorted(keep):f.write(f'{hashlib.sha256(q.read_bytes()).hexdigest()}  {q.relative_to(ROOT)}\n')
review=[q for q in keep if q.suffix!='.bin' and ('/runs/' not in str(q) or q.name in ['metrics.json','result.json','command.json'])]
with zipfile.ZipFile(OUT/'[OMEx]_metal_envelope_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for q in review:z.write(q,str(q.relative_to(ROOT)))
 z.write(OUT/'SHA256SUMS',str((OUT/'SHA256SUMS').relative_to(ROOT)))
with zipfile.ZipFile(OUT/'[OMEx]_metal_envelope_review.zip') as z:assert z.testzip() is None
print('Packaged',len(review),'review files; native histories remain in isolated runs/.')
