"""Static edge check only: no new edge trajectories."""
from pathlib import Path
import subprocess,csv,json
import numpy as np
from scipy.linalg import eigvals
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1];p=ROOT/'data/output/turn_on_mechanisms_20261117'
es=np.unique(np.r_[np.linspace(2.451,2.456,101),np.linspace(3.212,3.219,141)])
(p/'inputs/edge_energies.txt').write_text(''.join(f'{e:.17g}\n' for e in es))
subprocess.run([str(p/'turn_on_matrix.bin'),str(p/'tables/edge_matrices.tsv'),str(p/'inputs/edge_energies.txt'),'1.1'],check=True)
a=np.genfromtxt(p/'tables/edge_matrices.tsv',names=True,delimiter='\t');rows=[]
for r in a[a['N']==1]:
 A=np.array([r[f'Ar{i}']+1j*r[f'Ai{i}'] for i in range(9)]).reshape(3,3);v=eigvals(A);ix=v.real.argmax();rows.append(dict(energy_eV=r['energy_eV'],sigma_ps_inv=v[ix].real,imag_ps_inv=v[ix].imag))
with (p/'tables/edge_rates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
roots=[]
for lo,hi in [(2.451,2.456),(3.212,3.219)]:
 rr=[r for r in rows if lo<=r['energy_eV']<=hi];s=CubicSpline([r['energy_eV'] for r in rr],[r['sigma_ps_inv'] for r in rr]);r=s.roots(extrapolate=False);coarse=CubicSpline([x['energy_eV'] for x in rr[::2]],[x['sigma_ps_inv'] for x in rr[::2]]).roots(extrapolate=False);roots.append(dict(bracket_eV=[lo,hi],root_eV=float(r[0]),coarse_grid_root_eV=float(coarse[0]),grid_spacing_eV=float(rr[1]['energy_eV']-rr[0]['energy_eV'])));print(lo,hi,r)

(p/"tables/edge_roots.json").write_text(json.dumps(roots,indent=2))
