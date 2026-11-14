"""Reproduce the new autonomous maps and paired energies from frozen inputs."""
from pathlib import Path
import numpy as np,json,subprocess
p=Path(__file__).resolve().parent

def run(mode,out,*args):
 folder=p/out;folder.mkdir(parents=True,exist_ok=True);cmd=[str(p/'driver'),mode,str(folder),*map(str,args)];(folder/'command.json').write_text(json.dumps(cmd,indent=2))
 with (folder/'map.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
run('pair','inputs')
grid=np.loadtxt(p/'sources/Data/Pulling/output/ns_intensity_SS.dat')[:,0]
np.savetxt(p/'inputs/pulling_grid.txt',grid,fmt='%.17g');run('map','runs/pulling_map',1.2,p/'inputs/pulling_grid.txt')
G=json.loads((p/'inputs/grids.json').read_text())
for g,a in G.items():
 np.savetxt(p/f'inputs/map_grid_{g}.txt',a['frequencies'],fmt='%.17g');run('map',f'runs/frequency_map_{g}',g,p/f'inputs/map_grid_{g}.txt')
pair=json.loads((p/'inputs/pair.json').read_text());np.savetxt(p/'inputs/pulling_checks_grid.txt',[2.55,2.7,pair['threshold_eV'],2.95,3.1],fmt='%.17g');run('map','runs/pulling_check_map',1.2,p/'inputs/pulling_checks_grid.txt')
