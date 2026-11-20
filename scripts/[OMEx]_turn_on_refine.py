"""Small central-neighborhood runs; frozen corrected inputs and full-step intensities."""
from pathlib import Path
import subprocess,json,hashlib,gzip,shutil,concurrent.futures
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/turn_on_mechanisms_20261117';eth=2.8121971959422813
jobs=[]
for d in [0.,-.01,.01,-.001,.001,-.0001,.0001,-.00001,.00001,-.000001,.000001,2.8122-eth]:
 name='central_'+('rounded' if d==2.8122-eth else f'{d:+.0e}')
 jobs.append(dict(name=name,energy=eth+d,dt=1,T=80 if abs(d) in [0.,.001] else 40))
for d in [0.,-.001,.001,-.00001,.00001,-.000001,.000001,2.8122-eth]:jobs.append(dict(name='half_'+('rounded' if d==2.8122-eth else f'{d:+.0e}'),energy=eth+d,dt=.5,T=40))
(OUT/'inputs/targeted_jobs.json').write_text(json.dumps(jobs,indent=2))
def run(j):
 folder=OUT/'runs'/j['name'];folder.mkdir(exist_ok=True);record=folder/'result.json'
 config={**j,'gain':1.1,'E0_over_Esat':1e-9,'pump':'step','cadence_ps':.001,'initial_N':0,'initial_q':[0,0,0],'binary_sha256':hashlib.sha256((OUT/'gain_switching.bin').read_bytes()).hexdigest()}
 if record.exists():
  prior=json.loads(record.read_text());assert prior['config']==config and prior['returncode']==0;return
 cmd=[str(OUT/'gain_switching.bin'),'run',str(folder),'1.1',str(j['energy']),'1e-9','step',str(j['T']),str(j['dt']),'.001','1','0','1','0'];(folder/'command.json').write_text(json.dumps(cmd,indent=2))
 with (folder/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=f,timeout=120)
 record.write_text(json.dumps({'config':config,'returncode':r.returncode},indent=2));assert r.returncode==0
 with (folder/'intensity.f64').open('rb') as a,gzip.open(folder/'intensity.f64.gz','wb',compresslevel=1) as b:shutil.copyfileobj(a,b)
 (folder/'intensity.f64').unlink();print(j['name'],'complete',flush=True)
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(4) as ex:list(ex.map(run,jobs))
