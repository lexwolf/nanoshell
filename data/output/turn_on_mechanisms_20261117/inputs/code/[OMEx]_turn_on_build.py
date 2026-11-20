"""Build isolated diagnostics from the pinned library; never update the submodule."""
from pathlib import Path
import subprocess,json,hashlib,csv,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/turn_on_mechanisms_20261117'
for d in ['inputs','tables','runs','logs']:(OUT/d).mkdir(parents=True,exist_ok=True)
# Recover the same bounded diagnostic set from the original corrected campaign.
OLD=ROOT/'data/output/gain_switching_corrected_20261114';eth=2.8121971959422813
rows=[x for x in csv.DictReader((OLD/'tables/frequency_scan.csv').open()) if float(x['gain_ratio'])==1.1];selected=[]
for e in [2.48,2.50,2.55,2.6,2.7,2.798,2.8051,2.8122,2.8193,2.8264,2.9,3.,3.1,3.14,3.16]:
 x=min(rows,key=lambda x:abs(float(x['frequency_eV'])-e))
 if x not in selected:selected.append(x)
(OUT/'inputs/reused_selection.json').write_text(json.dumps(selected,indent=2))
energies=sorted(set([float(x['frequency_eV']) for x in selected]+[eth,2.5309774763480535,3.1201269537068566]+[eth+d for v in [1e-2,1e-3,1e-4,1e-5,1e-6] for d in [-v,v]]))
(OUT/'inputs/energies.txt').write_text(''.join(f'{e:.17g}\n' for e in energies))
for name in ['turn_on_matrix','gain_switching']:
 cmd=['g++','-O3','-std=c++17',*[f'-I{ROOT/"extern/nano_geo_matrix"/x}' for x in ['include','modules','modules/cup']],'-I/usr/include/eigen3',str(ROOT/f'src/[OMEx]_{name}.cxx'),'-o',str(OUT/f'{name}.bin'),'-lgsl','-lgslcblas','-lm','-larmadillo']
 (OUT/f'inputs/build_{name}.json').write_text(json.dumps(cmd,indent=2));subprocess.run(cmd,check=True)
prov={'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'library_HEAD':subprocess.check_output(['git','-C','extern/nano_geo_matrix','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'status':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),'compiler':subprocess.check_output(['g++','--version'],text=True),'sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [ROOT/'src/[OMEx]_gain_switching.cxx',ROOT/'src/[OMEx]_turn_on_matrix.cxx',OUT/'gain_switching.bin',OUT/'turn_on_matrix.bin']}}
(OUT/'inputs/build_provenance.json').write_text(json.dumps(prov,indent=2))
subprocess.run([str(OUT/'turn_on_matrix.bin'),str(OUT/'tables/matrices.tsv'),str(OUT/'inputs/energies.txt'),'1.1'],check=True)
