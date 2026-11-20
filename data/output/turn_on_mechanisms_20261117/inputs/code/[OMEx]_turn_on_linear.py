"""Tracked non-Hermitian modes and continuous-forcing seed diagnostics."""
from pathlib import Path
import itertools,json,csv
import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/output/turn_on_mechanisms_20261117'
def models():
 tab=np.genfromtxt(OUT/'tables/matrices.tsv',names=True,delimiter='\t');result={}
 for e in np.unique(tab['energy_eV']):
  a=tab[tab['energy_eV']==e];A=[];b=[]
  for r in a:
   A.append(np.array([r[f'Ar{i}']+1j*r[f'Ai{i}'] for i in range(9)]).reshape(3,3));b.append(np.array([r[f'br{i}']+1j*r[f'bi{i}'] for i in range(3)]))
  c=np.array([a[0][f'cr{i}']+1j*a[0][f'ci{i}'] for i in range(4)])
  result[float(e)]={'A0':A[0],'D':A[1]-A[0],'b0':b[0],'db':b[1]-b[0],'c':c,'tau1':float(a[0]['tau1_ps']),'Esat':float(a[0]['Esat']),'E0':float(a[0]['E0']),'k':float(a[0]['k_ps_inv'])}
 return result
def write(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def linear():
 summary=[];tracks=[];forcing=[]
 for e,m in models().items():
  A0,D,b0,db,c=[m[x] for x in ['A0','D','b0','db','c']];previous=None;track=[]
  for N in np.linspace(0,1,501):
   A=A0+N*D;val,V=eig(A)
   order=np.argsort(val.imag) if previous is None else min(itertools.permutations(range(3)),key=lambda p:np.sum(abs(val[list(p)]-previous)**2))
   val=val[list(order)];V=V[:,list(order)];previous=val.copy();W=np.linalg.inv(V);sep=min(abs(val[i]-val[j]) for i in range(3) for j in range(i))
   track.append(val)
   for j in range(3):tracks.append(dict(energy_eV=e,N=N,mode=j,sigma_ps_inv=val[j].real,imag_ps_inv=val[j].imag,condition=float(np.linalg.norm(V[:,j])*np.linalg.norm(W[j])),eigenvector_matrix_condition=float(np.linalg.cond(V)),min_eigenvalue_separation_ps_inv=sep))
  val=previous;idx=int(np.argmax(val.real));lam=val[idx];v=V[:,idx];left=W[idx];cv=c[:3]@v
  assert np.max(abs(W@V-np.eye(3)))<1e-10
  roots=[]
  for j in range(3):
   for k in range(500):
    if track[k][j].real*track[k+1][j].real<0:
     def fun(N):
      vs=eig(A0+N*D,right=False);return vs[np.argmin(abs(vs-track[k][j]))].real
     roots.append((brentq(fun,k/500,(k+1)/500),j))
  # Prescribed pump 1 but initial inversion 0: undepleted inversion ramp.
  def rhs(t,q):
   N=-np.expm1(-t/m['tau1']);return (A0+N*D)@q+b0+N*db
  sol=solve_ivp(rhs,[0,8],np.zeros(3,complex),method='DOP853',rtol=2e-10,atol=1e-24,dense_output=True)
  assert sol.success
  ts=np.linspace(0,8,4001);qs=sol.sol(ts);ds=c[:3]@qs+c[3]*m['E0'];Is=abs(ds/m['Esat'])**2
  # Include ALL non-Hermitian modal forcing and transient mixing; project after inversion settled.
  C=[]
  for t in [4.,6.,8.]:C.append(complex(cv*((left@sol.sol(t))+(left@(b0+db))/lam)*np.exp(-lam*t)))
  seed=abs(C[-1])/m['Esat'];frozen=abs(cv*(left@(b0+db))/lam)/m['Esat']
  frozenq=np.array([V@((np.exp(val*t)-1)*(W@(b0+db))/val) for t in ts]);frozenI=abs((frozenq@c[:3]+c[3]*m['E0'])/m['Esat'])**2
  np.savez_compressed(OUT/'tables'/f'linear_{e:.17g}.npz',t=ts,I=Is,q=qs,forced_N1_I=frozenI)
  summary.append(dict(energy_eV=e,mode=idx,sigma_ps_inv=lam.real,imag_ps_inv=lam.imag,instability_N=min([x[0] for x in roots],default=np.nan),mode_condition=float(np.linalg.norm(v)*np.linalg.norm(left)),V_condition=float(np.linalg.cond(V)),biorthogonal_residual=float(np.max(abs(W@V-np.eye(3)))),effective_output_seed_over_Esat=seed,frozen_N1_output_seed_over_Esat=frozen,seed_4_to_8_relative_change=abs(C[0]-C[-1])/abs(C[-1]),equivalent_ramp_delay_ps=np.log(frozen/seed)/lam.real,threshold_fixed_I_1e_4_prediction_ps=np.log(.01/seed)/lam.real,tau1_ps=m['tau1']))
 write('linear_summary.csv',summary);write('tracked_modes.csv',tracks)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':linear()
