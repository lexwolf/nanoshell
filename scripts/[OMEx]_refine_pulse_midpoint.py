"""Add a calculated gain=1.105 pulse case and verify it at half step."""
from pathlib import Path
import csv,gzip,hashlib,json,shutil,subprocess
import numpy as np
paths=__import__('[OMEx]_paths');p=paths.DATA
read=lambda q:dict(line.split('\t',1) for line in q.read_text().splitlines())
reports=[]
for name,dt in [('pulse_g1.105',1),('pulse_half_g1.105',.5)]:
 folder=p/'runs'/name;folder.mkdir(exist_ok=True)
 config=dict(name=name,gain=1.105,energy='threshold',seed=1e-9,pump='cos50',T=200,dt=dt,cadence=.002,factor=1,diagnostic=0,min_observe=0,source_sha256=hashlib.sha256(paths.DRIVER_SOURCE.read_bytes()).hexdigest())
 record=folder/'result.json'
 if record.exists():
  result=json.loads(record.read_text());assert result['configuration']==config and result['returncode']==0
 else:
  cmd=[str(paths.BINARY),'run',str(folder),'1.105','threshold','1e-9','cos50','200',str(dt),'.002','1','0','1','0']
  (folder/'command.json').write_text(json.dumps(cmd,indent=2))
  provenance={'repository_HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=paths.ROOT,text=True).strip(),'library_HEAD':subprocess.check_output(['git','-C','extern/nano_geo_matrix','rev-parse','HEAD'],cwd=paths.ROOT,text=True).strip(),'binary_sha256':hashlib.sha256(paths.BINARY.read_bytes()).hexdigest(),'worktree_status':subprocess.check_output(['git','status','--short'],cwd=paths.ROOT,text=True)}
  (folder/'provenance.json').write_text(json.dumps(provenance,indent=2))
  with (folder/'run.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=log,timeout=120)
  result=dict(configuration=config,returncode=r.returncode,parameters=read(folder/'parameters.tsv'),metrics=read(folder/'metrics.tsv'))
  record.write_text(json.dumps(result,indent=2));assert r.returncode==0
  with (folder/'intensity.f64').open('rb') as src,gzip.open(folder/'intensity.f64.gz','wb',compresslevel=1) as dst:shutil.copyfileobj(src,dst)
  (folder/'intensity.f64').unlink()
 cyc=np.atleast_1d(np.genfromtxt(folder/'cycles.tsv',names=True,delimiter='\t'));assert len(cyc)==4
 reports.append({k:float(cyc[-1][k]) for k in cyc.dtype.names})
a,b=reports
comparison={'gain_ratio':1.105,'baseline':a,'half_step':b,'peak_relative_difference':abs(a['peak_I']-b['peak_I'])/a['peak_I'],'peak_time_difference_ps':abs(a['peak_time_ps']-b['peak_time_ps']),'fwhm_difference_ps':abs(a['fwhm_ps']-b['fwhm_ps']),'metric_resolution':'every integration step; output cadence .002 ps'}
(p/'tables/pulse_midpoint_refinement.json').write_text(json.dumps(comparison,indent=2));print(json.dumps(comparison,indent=2))
result=json.loads((p/'runs/pulse_g1.105/result.json').read_text());params=result['parameters']
new=dict(gain_ratio=1.105,peak_I=a['peak_I'],peak_time_ps=a['peak_time_ps'],fwhm_ps=a['fwhm_ps'],half_rise_ps=a['half_rise_ps'],half_fall_ps=a['half_fall_ps'],Esat=float(params['Esat']),raw_E0=float(params['raw_E0']),status=result['metrics']['status'])
for filename,extra in [('gain_dependence.csv',[new]),('all_cycle_metrics.csv',[dict(gain_ratio=1.105,**{k:float(row[k]) for k in cyc.dtype.names}) for row in np.atleast_1d(np.genfromtxt(p/'runs/pulse_g1.105/cycles.tsv',names=True,delimiter='\t'))])]:
 q=p/'tables'/filename;rows=list(csv.DictReader(q.open()));fields=list(rows[0]);rows=[r for r in rows if float(r['gain_ratio'])!=1.105]+extra;rows.sort(key=lambda r:(float(r['gain_ratio']),float(r.get('cycle',0))))
 with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
rows=list(csv.DictReader((p/'tables/gain_dependence.csv').open()))
np.savetxt(p/'publication_plot_data/Pub_Ready/gain_dependence_linear_extended_plot.dat',[[float(r['gain_ratio']),float(r['fwhm_ps']),float(r['peak_I']),float(r['peak_time_ps']),0] for r in rows],fmt='%.17g',header='NEW all-step final-cycle metrics: gain FWHM peak_I/Isat peak_time unverified_flag')
__import__('[OMEx]_fit_gain_dependence').fit_gain_dependence()
subprocess.run(['gnuplot',str(paths.ROOT/'scripts/[OMEx]_gain_dependence_linear_extended.gp')],cwd='/tmp',check=True)
base=paths.IMAGES/'[OMEx]_gain_dependence_linear_extended'
subprocess.run(['pdftoppm','-png','-r','300','-singlefile',str(base)+'.pdf',str(base)],check=True)
