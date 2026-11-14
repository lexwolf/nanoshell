from pathlib import Path
import subprocess,concurrent.futures,json,time,hashlib,os,gzip,shutil,threading
p=Path(__file__).resolve().parent
G=json.loads((p/'inputs/grids.json').read_text());pair=json.loads((p/'inputs/pair.json').read_text())
jobs=[]
def add(name,g,e,seed,pump,T,dt=1,cad=.02,factor=1,diag=0,minobs=0):jobs.append(dict(name=name,gain=g,energy=e,seed=seed,pump=pump,T=T,dt=dt,cadence=cad,factor=factor,diagnostic=diag,min_observe=minobs))
for g,gr in G.items():
 for i,seed in enumerate(gr['probes']):add(f'probe_g{g}_{i:04d}',g,'threshold',seed,'step',100,cad=.02)
 old={}
 for row in gr['frequency_source_rows']:old[row[1]]=max(old.get(row[1],0),row[0])
 for i,e in enumerate(gr['frequencies']):
  t=old[e];T=min(2000,max(240,1.5*t+180));add(f'freq_g{g}_{i:04d}',g,e,1e-9,'step',T,diag=1,minobs=1.2*t+40)
gains=sorted(set([round(1.05+i*.01,2) for i in range(25)]+[1.3,1.33,1.35,1.38,1.4,1.43,1.45,1.48,1.5]))
for g in gains:add(f'pulse_g{g:.2f}',g,'threshold',1e-9,'cos50',200,cad=.002)
for side in ['left','right']:
 for dt in [1,.5]:add(f'pair_{side}_dt{dt:g}',1.1,pair[side+'_eV'],1e-9,'step',200,dt=dt,cad=5*6.856369348370542e-5)
add('pair_left_dense',1.1,pair['left_eV'],1e-9,'step',200,cad=6.856369348370542e-5)
add('pair_right_literal_1p11',1.1,1.11*pair['threshold_eV'],1e-9,'step',200,cad=.001)
for g in [1.05,1.1,1.15,1.2,1.25,1.5]:add(f'pulse_half_g{g:.2f}',g,'threshold',1e-9,'cos50',200,dt=.5,cad=.002)
for e in [2.55,2.7,'threshold',2.95,3.1]:add(f'pulling_check_{e}',1.2,e,1e-8,'step',200,cad=.002)
(p/'inputs/jobs.json').write_text(json.dumps(jobs,indent=2));sha=hashlib.sha256((p/'driver.cxx').read_bytes()).hexdigest();lock=threading.Lock();done=0;t0=time.monotonic()
def read(q):return dict(line.split('\t',1) for line in q.read_text().splitlines())
def run(j):
 global done
 folder=p/'runs'/j['name'];folder.mkdir(exist_ok=True);record=folder/'result.json'
 config={**j,'source_sha256':sha}
 if record.exists():
  prior=json.loads(record.read_text())
  if prior.get('configuration')!=config:raise RuntimeError(f'Configuration mismatch at {folder}; use a new isolated campaign directory')
  if prior.get('returncode')==0:return prior
 args=[j['gain'],j['energy'],j['seed'],j['pump'],j['T'],j['dt'],j['cadence'],j['factor'],j['diagnostic'],1,j['min_observe']]
 cmd=[str(p/'driver'),'run',str(folder),*map(str,args)];(folder/'command.json').write_text(json.dumps(cmd,indent=2))
 with (folder/'run.log').open('w') as f:
  try:r=subprocess.run(cmd,stdout=f,stderr=f,timeout=120);rc=r.returncode
  except subprocess.TimeoutExpired:rc=-999
 row={'configuration':config,'returncode':rc,'parameters':read(folder/'parameters.tsv') if (folder/'parameters.tsv').exists() else {},'metrics':read(folder/'metrics.tsv') if (folder/'metrics.tsv').exists() else {}}
 if (folder/'intensity.f64').exists():
  with (folder/'intensity.f64').open('rb') as src,gzip.open(folder/'intensity.f64.gz','wb',compresslevel=1) as dst:shutil.copyfileobj(src,dst)
  (folder/'intensity.f64').unlink()
 record.write_text(json.dumps(row,indent=2))
 with lock:
  done+=1
  if done%100==0:print(f'{done}/{len(jobs)} completed; elapsed {time.monotonic()-t0:.1f}s',flush=True)
 return row
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(int(os.environ.get('REGEN_WORKERS','16'))) as ex:results=list(ex.map(run,jobs))
 (p/'tables/all_results.json').write_text(json.dumps(results,indent=2));print('BATCH_FINISHED',len(results),time.monotonic()-t0,flush=True)
