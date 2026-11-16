from pathlib import Path
import shutil,json,subprocess,hashlib,numpy as np
r=Path('/home/alessandro/GitHub/Academia/GEO_MATRIX');o=r/'nanoshell/data/output/gain_switching_corrected_20261114';h=Path('/home/alessandro/Nextcloud/Work/Research/Papers/Active/Gain Switching')
for sub in ['gnuplot','scripts']:
 shutil.copytree(h/'work/Gain_Switching'/sub,o/'sources'/sub,dirs_exist_ok=True)
shutil.copytree(h/'round 1/LaTeX',o/'sources/LaTeX',dirs_exist_ok=True)
shutil.copytree(h/'round 1/docs/residual_modulation_diagnostic',o/'sources/residual_modulation_diagnostic',dirs_exist_ok=True)
shutil.copy2('/home/alessandro/Documents/ChatGPT/Plasmonic Nano-Shell Gain Switching/normalization_figure_audit.md',o/'sources/normalization_figure_audit.md')
for sub in ['ProbeSweep/middle','FreqSweep','Pulling']:
 shutil.copytree(h/'work/Gain_Switching/Data'/sub,o/'sources/Data'/sub,dirs_exist_ok=True)
for repo in ['nanoshell','nano_geo_matrix','nanoshell/extern/nano_geo_matrix']:
 cmd=['git','-C',str(r/repo)];(o/'sources'/(repo.replace('/','_')+'_state.txt')).write_text(''.join(subprocess.check_output(cmd+x,text=True) for x in [['show','-s','--format=fuller','HEAD'],['status','--porcelain'],['remote','-v']]))
shutil.copy2(r/'nanoshell/src/residual_diagnostic.hpp',o/'residual_diagnostic.hpp')
grids={}
for label,g in [('1.05',1.05),('1.1',1.1),('1.15',1.15)]:
 a=np.loadtxt(o/f'sources/Data/ProbeSweep/middle/{label}Gth.dat');b=np.loadtxt(o/f'sources/Data/FreqSweep/{label}Gth/{label}Gth.dat')
 grids[str(g)]={'probes':a[:,2].tolist(),'frequencies':sorted(set(b[:,1])),'frequency_source_rows':b.tolist()}
(o/'inputs/grids.json').write_text(json.dumps(grids,indent=2))
(o/'inputs/parameters.json').write_text(json.dumps({'a_nm':10,'rho':.6,'Dome_eV':.15,'metal':'silver','model':'drude','active':'lorentz','core':'silica','exterior':'water','initial_q':[0,0,0],'initial_N':0,'gain_center':'computed_threshold','pump_step':'Ntarget=1 for t>=0','pump_periodic':'Ntarget=(1-cos(2*pi*t/50))/2','tau1_policy':'require_match','scan_max_duration_ps':2000,'probe_max_duration_ps':100,'gain_cycles_duration_ps':200,'normalization_gain':'fixed_Gmax'} ,indent=2))
