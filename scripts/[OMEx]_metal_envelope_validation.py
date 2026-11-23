"""Separate exact-metal-envelope validation; no changes to the legacy solver."""
from pathlib import Path
import json,csv,hashlib,subprocess,gzip,shutil,concurrent.futures,sys
import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import brentq,root
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/metal_envelope_validation_20261120';IMG=ROOT/'img/output/[OMEx]_metal_envelope_validation'
def constants():return json.loads((OUT/'inputs/constants.json').read_text())
def write(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def matrix(N,G,wg,wr):
 c=constants();p=np.array(c['p']);tau2=c['tau2'];gam=c['gamma_normalized'];A=np.zeros((5,5),complex);b=np.zeros(5,complex)
 A[0,:3]=-1j*G/tau2*N*p[0,:3];A[0,0]+=1j*(wr-wg)-1/tau2;b[0]=-1j*G/tau2*N*p[0,3]
 for j in [1,2]:A[j,j]=1j*wr;A[j,j+2]=1;A[j+2,:3]=p[j,:3];A[j+2,j+2]=1j*wr-2*gam;b[j+2]=p[j,3]
 return A,b

def build():
 for d in ['inputs','runs','tables','logs']:(OUT/d).mkdir(parents=True,exist_ok=True)
 cmd=['g++','-O3','-std=c++17',*[f'-I{ROOT/"extern/nano_geo_matrix"/x}' for x in ['include','modules','modules/cup']],'-I/usr/include/eigen3',str(ROOT/'src/[OMEx]_metal_envelope_comparison.cxx'),'-o',str(OUT/'comparison.bin'),'-lgsl','-lgslcblas','-lm','-larmadillo']
 (OUT/'inputs/build.json').write_text(json.dumps(cmd,indent=2));subprocess.run(cmd,check=True);subprocess.run([str(OUT/'comparison.bin'),'constants',str(OUT/'inputs/constants.json')],check=True)

def spectral():
 import sympy as sp
 t,w,g=sp.symbols('t w gamma',real=True);P=sp.Function('P')(t)
 result=sp.expand((sp.diff(P*sp.exp(-sp.I*w*t),t,2)+2*g*sp.diff(P*sp.exp(-sp.I*w*t),t))*sp.exp(sp.I*w*t))
 target=sp.diff(P,t,2)+2*(g-sp.I*w)*sp.diff(P,t)-(w*w+2*sp.I*g*w)*P;assert sp.simplify(result-target)==0
 (OUT/'inputs/envelope_sign_check.txt').write_text(str(result)+'\nSymbolic residual: '+str(sp.simplify(result-target)))
 c=constants();ep=c['plasma_eV'];eth=c['old_threshold_eV'];wg=eth/ep;G=1.1*c['old_Gth'];p=np.array(c['p']);k=c['k_ps_inv'];wr=wg
 # Same physical G and gain center; dominant autonomous eigenmode of the unforced problem.
 def dom(N,G=G,center=wg,ref=wr):
  v=eigvals(matrix(N,G,center,ref)[0]);return v[np.argmax(v.real)]
 Gcrit=brentq(lambda g:dom(1,g).real,.05,.3,xtol=1e-14)
 aligned=root(lambda x:[dom(1,x[0],x[1]/ep,x[1]/ep).real,dom(1,x[0],x[1]/ep,x[1]/ep).imag],[Gcrit,eth],tol=1e-11);assert aligned.success
 fixed=dom(1,Gcrit);sat=dom(Gcrit/G)
 summary=dict(old_Gth=c['old_Gth'],old_threshold_eV=eth,fixed_gain_center_Gth=Gcrit,fixed_gain_center_threshold_emission_eV=eth-ep*fixed.imag,aligned_Gth=float(aligned.x[0]),aligned_gain_center_eV=float(aligned.x[1]),G=G,full_inversion_sigma_ps_inv=float(dom(1).real*k),full_inversion_nu_ps_inv=float(dom(1).imag*k),critical_N=Gcrit/G,saturated_emission_eV=eth-ep*sat.imag)
 # Nonlinear stationary intensity from the null eigenmode and the unchanged inversion equation.
 from scipy.linalg import eig
 val,V=eig(matrix(Gcrit/G,G,wg,wr)[0]);j=np.argmin(abs(val.real));v=V[:,j];den=np.imag(v[0]*np.conj(p[0,:3]@v[:3]));a2=-(1-Gcrit/G)/(c['tau1']*den)
 # Esat will be read from the existing shared-utility snapshot for exactly this G.
 old=ROOT/'data/output/gain_switching_corrected_20261114/runs/pair_left_dt1/parameters.tsv';par=dict(x.split('\t',1) for x in old.read_text().splitlines());Isat=float(par['Isat']);summary['stationary_J_at_identical_G']=float(a2*abs(p[3,:3]@v[:3])**2/Isat)
 (OUT/'tables/thresholds.json').write_text(json.dumps(summary,indent=2))
 passive=[];cov=[]
 for e in [2.53,eth,3.12]:
  wpr=e/ep;gamma=c['gamma_normalized'];epsm=c['eps_inf']-1/(wpr*(wpr+2j*gamma));eb=c['eps_core'];ee=c['eps_host'];rho3=.6**3
  analytic=((epsm-ee)*(eb+2*epsm)+rho3*(eb-epsm)*(ee+2*epsm))/((epsm+2*ee)*(eb+2*epsm)+2*rho3*(epsm-ee)*(eb-epsm))
  for ref in [eth/ep,(eth+.1)/ep]:
   A,b=matrix(0,G,wg,ref);delta=1j*(ref-wpr);q=np.linalg.solve(delta*np.eye(5)-A,b);dip=p[3,:3]@q[:3]+p[3,3];chi=-1/(wpr*(wpr+2j*gamma));err=max(abs(q[i]-chi*(p[i,:3]@q[:3]+p[i,3])) for i in [1,2]);passive.append(dict(probe_eV=e,reference_eV=ref*ep,analytic_Re=analytic.real,analytic_Im=analytic.imag,matrix_Re=dip.real,matrix_Im=dip.imag,relative_dipole_error=abs(dip-analytic)/abs(analytic),metal_constitutive_abs_error=err))
 for ref in [eth-.3,eth,eth+.1,eth+.3]:
  A,_=matrix(1,G,wg,ref/ep);A0,_=matrix(1,G,wg,wg);v=eigvals(A);v0=eigvals(A0);shift=1j*(ref-eth)/ep;res=np.linalg.norm(A-A0-shift*np.eye(5));cov.append(dict(reference_eV=ref,matrix_covariance_residual=res,max_growth_ps_inv=max(v.real)*k,growth_difference_ps_inv=(max(v.real)-max(v0.real))*k))
 write('passive_frequency_validation.csv',passive);write('carrier_spectral_validation.csv',cov)
 assert max(x['relative_dipole_error'] for x in passive)<1e-11
 assert max(abs(x['growth_difference_ps_inv']) for x in cov)<1e-9
 print(json.dumps(summary,indent=2))

def runjobs(jobs):
 c=constants();ep=c['plasma_eV'];eth=c['old_threshold_eV'];G=1.1*c['old_Gth'];binary=OUT/'comparison.bin';hashbin=hashlib.sha256(binary.read_bytes()).hexdigest()
 def run(j):
  j={'G':G,'center':eth,'reference':j.get('probe',eth),'T':40,'dt':.125,'prep':'rest','cadence':.001,**j};f=OUT/'runs'/j['name'];f.mkdir(exist_ok=True);cfg={**j,'binary_sha256':hashbin};record=f/'result.json'
  if record.exists():
   r=json.loads(record.read_text());assert r['config']==cfg and r['returncode']==0;return
  cmd=[str(binary),str(f),*[str(j[k]) for k in ['G','center','probe','reference','T','dt','prep','cadence']]];(f/'command.json').write_text(json.dumps(cmd,indent=2))
  with (f/'run.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=log,timeout=180)
  record.write_text(json.dumps(dict(config=cfg,returncode=r.returncode),indent=2));assert r.returncode==0
  with (f/'intensity.f64').open('rb') as a,gzip.open(f/'intensity.f64.gz','wb',compresslevel=1) as b:shutil.copyfileobj(a,b)
  (f/'intensity.f64').unlink();print(j['name'],'completed',flush=True)
 with concurrent.futures.ThreadPoolExecutor(3) as ex:list(ex.map(run,jobs))

def pilot():
 c=constants();eth=c['old_threshold_eV'];jobs=[]
 for dt in [1,.5,.25,.125,.0625]:
  for name,e in [('threshold',eth),('left',2.5309774763480535),('right',3.1201269537068566)]:
   if dt in [1,.5] and name!='threshold':continue
   jobs.append(dict(name=f'{name}_dt{dt:g}',probe=e,dt=dt,T=30))
 for dt in [.125,.0625]:jobs.append(dict(name=f'carrier_dt{dt:g}',probe=eth,reference=eth+.1,dt=dt,T=30))
 for e,name in [(2.53,'left'),(eth,'center'),(3.12,'right')]:jobs.append(dict(name=f'passive_{name}',probe=e,reference=eth+.1,prep='passive',dt=.0625,T=2,cadence=.0001))
 runjobs(jobs)

def central():
 assert json.loads((OUT/'tables/pilot_validation.json').read_text())['passed']
 c=constants();eth=c['old_threshold_eV'];th=json.loads((OUT/'tables/thresholds.json').read_text());jobs=[]
 for d in [-.01,-.001,-.0001,-.00001,.00001,.0001,.001,.01]:jobs.append(dict(name=f'central_{d:+.0e}',probe=eth+d,T=30,dt=.125))
 for d in [-.001,.001]:jobs.append(dict(name=f'centralhalf_{d:+.0e}',probe=eth+d,T=30,dt=.0625))
 for name,e in [('threshold',eth),('left',2.5309774763480535),('right',3.1201269537068566)]:jobs.append(dict(name='slaved_'+name,probe=e,prep='slaved',dt=.125,T=30))
 # The full autonomous map is carrier invariant; use an explicitly specified symmetric physical-probe pair, not an assumed old folding relation.
 newcenter=th['aligned_gain_center_eV'];left=.9*newcenter;right=2*th['saturated_emission_eV']-left
 for name,e in [('threshold',newcenter),('left',left),('right',right)]:jobs.append(dict(name='matched_'+name,probe=e,G=1.1*th['aligned_Gth'],center=newcenter,dt=.125,T=30))
 (OUT/'inputs/matched_pair_definition.json').write_text(json.dumps(dict(left_eV=left,right_eV=right,common_autonomous_emission_eV=th['saturated_emission_eV'],definition='Symmetric physical probe detunings about the computed common autonomous frequency. Not a unique two-branch frequency-pulling pair: the restored autonomous map is constant in carrier reference.',G=1.1*th['aligned_Gth'],gain_center_eV=newcenter),indent=2))
 # Explicit threshold-aligned comparison if the computed threshold changes materially.
 (OUT/'inputs/targeted_central_jobs.json').write_text(json.dumps(jobs,indent=2));runjobs(jobs)

if __name__=='__main__':globals()[sys.argv[1]]()
