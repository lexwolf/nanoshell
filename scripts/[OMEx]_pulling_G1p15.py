"""Isolated Fig. 1 candidate at gain 1.15; never overwrite approved assets."""
from pathlib import Path
import subprocess,json,hashlib,gzip,shutil,csv,concurrent.futures
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/pulling_G1p15_20261120';OLD=ROOT/'data/output/gain_switching_corrected_20261114';IMG=ROOT/'img/output/[OMEx]_gain_switching';SCRIPTS=ROOT/'scripts'
for sub in ['inputs','runs','tables','plot_data','logs']:(OUT/sub).mkdir(parents=True,exist_ok=True)
read=lambda q:dict(line.split('\t',1) for line in q.read_text().splitlines())
def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
approved=[IMG/'[OMEx]_pulling.pdf',IMG/'[OMEx]_pulling.png',SCRIPTS/'[OMEx]_pulling.gp']
before={str(q.relative_to(ROOT)):sha(q) for q in approved}
source=ROOT/'src/[OMEx]_gain_switching.cxx';binary=OUT/'driver.bin'
cmd=['g++','-O3','-std=c++17',*[f'-I{ROOT/"extern/nano_geo_matrix"/x}' for x in ['include','modules','modules/cup']],'-I/usr/include/eigen3',str(source),'-o',str(binary),'-lgsl','-lgslcblas','-lm','-larmadillo']
(OUT/'inputs/build.json').write_text(json.dumps(cmd,indent=2));subprocess.run(cmd,check=True)
prov={'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'library_HEAD':subprocess.check_output(['git','-C','extern/nano_geo_matrix','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'status_before':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),'compiler':subprocess.check_output(['g++','--version'],text=True),'source_sha256':sha(source),'binary_sha256':sha(binary),'approved_sha256_before':before}
(OUT/'inputs/provenance.json').write_text(json.dumps(prov,indent=2));shutil.copy2(source,OUT/'inputs/driver.cxx')
eth=json.loads((OLD/'inputs/pair.json').read_text())['threshold_eV'];grid=np.loadtxt(OLD/'inputs/pulling_grid.txt');grid=np.unique(np.r_[grid,eth]);np.savetxt(OUT/'inputs/grid.txt',grid,fmt='%.17g')
energies=[2.55,2.7,eth,2.95,3.1];np.savetxt(OUT/'inputs/check_energies.txt',energies,fmt='%.17g')
def maprun(name,grid):
 folder=OUT/'runs'/name;folder.mkdir(exist_ok=True);cmd=[str(binary),'map',str(folder),'1.15',str(grid)];(folder/'command.json').write_text(json.dumps(cmd,indent=2))
 with (folder/'run.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
maprun('stationary',OUT/'inputs/grid.txt');maprun('check_map',OUT/'inputs/check_energies.txt')
def run(pair):
 i,e=pair;folder=OUT/'runs'/f'check_{i}';folder.mkdir(exist_ok=True)
 cmd=[str(binary),'run',str(folder),'1.15',str(e),'1e-8','step','200','1','.002','1','0','1','0'];(folder/'command.json').write_text(json.dumps(cmd,indent=2))
 with (folder/'run.log').open('w') as log:res=subprocess.run(cmd,stdout=log,stderr=log,timeout=120)
 (folder/'result.json').write_text(json.dumps({'returncode':res.returncode},indent=2));assert res.returncode==0
 if (folder/'intensity.f64').exists():
  with (folder/'intensity.f64').open('rb') as a,gzip.open(folder/'intensity.f64.gz','wb',compresslevel=1) as b:shutil.copyfileobj(a,b)
  (folder/'intensity.f64').unlink()
 print('Time-domain check complete:',e,flush=True)
with concurrent.futures.ThreadPoolExecutor(4) as ex:list(ex.map(run,enumerate(energies)))
pars=read(OUT/'runs/check_0/parameters.tsv');Isat=float(pars['Isat']);Esat=float(pars['Esat'])
assert abs(float(pars['E0_over_Esat'])/1e-8-1)<1e-12 and pars['supplied_tau1_matches']=='1'
a=np.genfromtxt(OUT/'runs/stationary/map.tsv',names=True,delimiter='\t');valid=(a['I_raw']>0)&(abs(a['Omega_normalized'])>0);positive=a[a['I_raw']>0];peak=float(max(a['I_raw'])/Isat)
np.savetxt(OUT/'plot_data/ns_intensity_SS.dat',np.c_[a['energy_eV'],a['Omega_normalized']*9.6,np.zeros(len(a)),np.zeros(len(a)),a['I_raw'],a['Nss'],np.zeros(len(a))],fmt='%.17g',header='NEW G/Gth=1.15 autonomous stationary construction; energy Omega_eV unused unused raw_I Nss unused')
em=a['emitted_eV'][valid];I=a['I_raw'][valid];ix=np.argsort(em);shifted=np.vstack([[1.8,0],[em.min(),0],np.c_[em[ix],I[ix]],[em.max(),0],[3.6,0]])
np.savetxt(OUT/'plot_data/ome_p3.dat',shifted,fmt='%.17g',header='NEW G/Gth=1.15 frequency-pulled stationary raw intensity; artificial closure points outside the mapped band')
(OUT/'plot_data/omeB.dat').write_text(f'{eth:.17g}\n');(OUT/'plot_data/nanosphere_eV.dat').write_text(f'10 .15 {eth:.17g} {pars["Gmax"]} 1.8 3.6 silver drude lorentz water {pars["raw_E0"]} .6 silica\n')
levels=np.array([.17,.566,.884])*peak*100
(OUT/'plot_data/parameters.gp').write_text(f'candidate_Isat = {Isat:.17g}\ncandidate_Esat = {Esat:.17g}\n')
# Copy approved styling, replacing only candidate data/outputs and gain-dependent arrow levels/range.
s=(SCRIPTS/'[OMEx]_pulling.gp').read_text();s=s.replace("load 'scripts/[OMEx]_publication_parameters.gp'", "load 'data/output/pulling_G1p15_20261120/plot_data/parameters.gp'")
s=s.replace('omex_Isat_120','candidate_Isat')
base='data/output/gain_switching_corrected_20261114/publication_plot_data/Data/Pulling'
for old,new in [(f'{base}/output/ns_intensity_SS.dat','ns_intensity_SS.dat'),(f'{base}/output/pulling/ome_p3.dat','ome_p3.dat'),(f'{base}/output/pulling/omeB.dat','omeB.dat'),(f'{base}/input/nanosphere_eV.dat','nanosphere_eV.dat')]:s=s.replace(old,f'data/output/pulling_G1p15_20261120/plot_data/{new}')
s=s.replace('array arrow_levels[3] = [1.44, 4.80, 7.488]','array arrow_levels[3] = ['+', '.join(f'{x:.17g}' for x in levels)+']')
s=s.replace('set yrange [0:10]','set yrange [0:7]')
s=s.replace('[OMEx]_pulling.pdf','[OMEx]_pulling_G1p15.pdf');q=SCRIPTS/'[OMEx]_pulling_G1p15.gp';q.write_text(s)
with (OUT/'logs/gnuplot.log').open('w') as log:subprocess.run(['gnuplot',str(q)],cwd='/tmp',check=True,stdout=log,stderr=log)
b=IMG/'[OMEx]_pulling_G1p15';subprocess.run(['pdftoppm','-png','-r','300','-singlefile',str(b)+'.pdf',str(b)],check=True)
checks=[];prediction=np.genfromtxt(OUT/'runs/check_map/map.tsv',names=True,delimiter='\t')
for i,e in enumerate(energies):
 folder=OUT/'runs'/f'check_{i}';tr=np.loadtxt(folder/'trajectory.tsv');tail=tr[tr[:,0]>=tr[-1,0]-20];r=prediction[i];pred=float(r['I_raw']/Isat);meas=float(tail[:,1].mean());pars_i=read(folder/'parameters.tsv');metrics=read(folder/'metrics.tsv')
 checks.append(dict(energy_eV=e,stationary_I_over_Isat=pred,tail20ps_mean_I_over_Isat=meas,relative_mean_difference=meas/pred-1,tail20ps_peak_to_peak_over_mean=float(np.ptp(tail[:,1])/meas),stationary_emitted_eV=float(r['emitted_eV']),measured_phase_emitted_eV=float(metrics['phase_emitted_eV']),stationary_N=float(r['Nss']),tail20ps_mean_N=float(tail[:,4].mean()),status=metrics['status'],raw_E0=float(pars_i['raw_E0']),E0_over_Esat=float(pars_i['E0_over_Esat'])))
with (OUT/'tables/time_checks.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
summary={'gain_ratio':1.15,'parameters':pars,'stationary_grid_samples':len(a),'positive_excitation_sample_bounds_eV':[float(positive['energy_eV'].min()),float(positive['energy_eV'].max())],'emission_sample_bounds_eV':[float(em.min()),float(em.max())],'maximum_I_over_Isat':peak,'arrow_levels_display':levels.tolist(),'checks':checks,'approved_unchanged':all(sha(ROOT/path)==hash for path,hash in before.items())}
assert summary['approved_unchanged'];(OUT/'tables/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

# Record the interpolation accuracy of the newly calculated paired arrows.
arrow_rows=[]
for level in levels:
 hits=[]
 for i in range(1,len(a)):
  y0=a['I_raw'][i-1]/Isat*100;y1=a['I_raw'][i]/Isat*100
  if (y0-level)*(y1-level)<0:
   w=(level-y0)/(y1-y0);hits.append((float(a['energy_eV'][i-1]+w*(a['energy_eV'][i]-a['energy_eV'][i-1])),float(a['emitted_eV'][i-1]+w*(a['emitted_eV'][i]-a['emitted_eV'][i-1]))))
 assert len(hits)==2
 arrow_rows.append(dict(level_I_over_Isat=float(level/100),left_energy_eV=hits[0][0],right_energy_eV=hits[1][0],left_mapped_emission_eV=hits[0][1],right_mapped_emission_eV=hits[1][1],interpolation_pair_mismatch_eV=hits[1][1]-hits[0][1]))
with (OUT/'tables/arrow_coordinates.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(arrow_rows[0]),lineterminator='\n');w.writeheader();w.writerows(arrow_rows)
# Snapshot the exact source and candidate plotting script used in this run.
for src in [SCRIPTS/'[OMEx]_pulling_G1p15.py',q,ROOT/'src/residual_diagnostic.hpp']:shutil.copy2(src,OUT/'inputs'/src.name)
