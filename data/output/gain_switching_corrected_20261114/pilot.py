from pathlib import Path
import subprocess,concurrent.futures,json,time,numpy as np
p=Path(__file__).resolve().parent
jobs=[('constant_ref',1.1,'threshold',1e-9,'step',100,1,.002,1,0,0),('constant_fast',1.1,'threshold',1e-9,'step',100,1,.002,1,0,1),('constant_half',1.1,'threshold',1e-9,'step',100,.5,.002,1,0,1),('periodic_ref',1.1,'threshold',1e-9,'cos50',200,1,.002,1,0,0),('periodic_fast',1.1,'threshold',1e-9,'cos50',200,1,.002,1,0,1),('periodic_half',1.1,'threshold',1e-9,'cos50',200,.5,.002,1,0,1),('edge',1.05,2.51003,1e-9,'step',2000,1,.05,1,1,1),('strong',1.1,'threshold',1.,'step',100,1,.002,1,1,1),('strong_half',1.1,'threshold',1.,'step',100,.5,.002,1,0,1)]
def run(j):
 name,*args=j;folder=p/'pilot'/name
 cmd=[str(p/'driver'),'run',str(folder),*map(str,args),"0"];folder.mkdir(exist_ok=True);(folder/'command.json').write_text(json.dumps(cmd));t=time.monotonic()
 with (folder/'log.txt').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=log)
 print(name,r.returncode,time.monotonic()-t,flush=True)
with concurrent.futures.ThreadPoolExecutor(8) as ex:list(ex.map(run,jobs))
report={}
for k in ['constant','periodic']:
 a=np.fromfile(p/'pilot'/f'{k}_ref/intensity.f64');b=np.fromfile(p/'pilot'/f'{k}_fast/intensity.f64');c=np.fromfile(p/'pilot'/f'{k}_half/intensity.f64')
 report[k]={'reference_fast_bit_identical':bool(np.array_equal(a,b)),'refinement_max_I_difference':float(np.max(np.abs(a-c[::2][:len(a)]))),'refinement_relative_to_peak':float(np.max(np.abs(a-c[::2][:len(a)]))/max(a))}
(p/'pilot/results.json').write_text(json.dumps(report,indent=2));print(report)
