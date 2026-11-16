"""Export corrected numerical results to the publication gnuplot layouts; no simulations."""
from pathlib import Path
import csv,json,hashlib,subprocess,shutil,numpy as np
paths=__import__('[OMEx]_paths');p=paths.DATA;root=paths.ROOT;scripts=root/'scripts';out=p/'publication_plot_data'
paths.IMAGES.mkdir(parents=True,exist_ok=True)
def rows(name):return list(csv.DictReader((p/'tables'/name).open()))
def folder(n):return p/'runs'/n
def record(n):return json.loads((folder(n)/'result.json').read_text())
def write(rel,a,header):
 q=out/rel;q.parent.mkdir(parents=True,exist_ok=True);np.savetxt(q,a,fmt='%.17g',header=header)
def trajectory(n,rel):
 a=np.loadtxt(folder(n)/'trajectory.tsv');r=record(n);E=float(r['parameters']['raw_E0']);dip=a[:,2]+1j*a[:,3]
 assert np.allclose(np.abs(dip)**2/float(r['parameters']['Isat']),a[:,1],rtol=1e-12,atol=0), n
 write(rel,np.c_[a[:,0],a[:,2]/E,a[:,3]/E,a[:,2],a[:,3],np.abs(dip),a[:,4]],f'NEW corrected run {n}; energy_eV={r["parameters"]["frequency_eV"]}; Esat={r["parameters"]["Esat"]}; E0/Esat={r["parameters"]["E0_over_Esat"]}\ncolumns: t Re_alpha Im_alpha Re_dip Im_dip modulus N')
 return a
probe=rows('probe_scan.csv');freq=rows('frequency_scan.csv');pulse=rows('gain_dependence.csv');cycles=rows('all_cycle_metrics.csv');pair=rows('paired_comparison.csv');pair_config=json.loads((p/'inputs/pair.json').read_text())
scales={g:float(record(f'pulse_g{g:.2f}')['parameters']['Esat']) for g in [1.05,1.1,1.15]}
for g,token in [(1.05,'1.05Gth'),(1.1,'1.1Gth'),(1.15,'1.15Gth')]:
 r=[x for x in probe if float(x['gain_ratio'])==g];e=scales[g]
 write(f'Data/ProbeSweep/middle/{token}.dat',[[float(x['t_on_ps']),float(x['raw_E0']),float(x['E0_over_Esat']),np.sqrt(float(x['maximum_I']))*e] for x in r],f'NEW corrected probes; Esat={e}; t E0 E0/Esat dipmax; all-step metrics')
 r=[x for x in freq if float(x['gain_ratio'])==g];a=[]
 for x in r:
  energy=float(x['frequency_eV']);em=float(x['emitted_phase_eV']);ton=float(x['t_on_ps']) if x['t_on_ps'] else np.nan
  a.append([ton,energy,energy/9.6,(energy-em)/9.6,em,np.sqrt(float(x['maximum_I']))*e])
 a=np.array(a);header='NEW corrected frequency scan; t_on or nan, excitation_eV, excitation_normalized, phase_shift_normalized, measured_emission_eV, dipmax; incomplete timings are nan'
 write(f'Data/FreqSweep/{token}/{token}.dat',a,header)
 write(f'Data/FreqSweep/{token}/Left.dat',a[a[:,1]<pair_config['threshold_eV']],header)
 write(f'Data/FreqSweep/{token}/Right.dat',a[a[:,1]>=pair_config['threshold_eV']],header)
 trajectory(f'pulse_g{g:.2f}',f'Data/GainSwitching/{token}/omegath.dat')
trajectory('pair_left_dt1','Data/TimeBehavior/0.9omegath.dat');trajectory('pair_right_dt1','Data/TimeBehavior/1.11omegath.dat')
# Fresh autonomous construction, not rescaled historical ordinates.
a=np.genfromtxt(folder('pulling_map')/'map.tsv',names=True,delimiter='\t');r=record('pulling_check_threshold');Isat=float(r['parameters']['Isat']);valid=(a['I_raw']>0)&(abs(a['Omega_normalized'])>0)
write('Data/Pulling/output/ns_intensity_SS.dat',np.c_[a['energy_eV'],a['Omega_normalized']*9.6,np.zeros(len(a)),np.zeros(len(a)),a['I_raw'],a['Nss'],np.zeros(len(a))],'NEW autonomous map: energy Omega_eV unused unused raw_I N unused')
em=a['emitted_eV'][valid];I=a['I_raw'][valid];ix=np.argsort(em);shifted=np.vstack([[1.8,0],[em.min(),0],np.c_[em[ix],I[ix]],[em.max(),0],[3.6,0]])
write('Data/Pulling/output/pulling/ome_p3.dat',shifted,'NEW autonomous frequency-pulled raw intensity')
q=out/'Data/Pulling/output/pulling/omeB.dat';q.write_text(f'{pair_config["threshold_eV"]:.17g}\n')
q=out/'Data/Pulling/input/nanosphere_eV.dat';q.parent.mkdir(parents=True,exist_ok=True);q.write_text(f'10 .15 {pair_config["threshold_eV"]:.17g} {float(r["parameters"]["Gmax"]):.17g} 1.8 3.6 silver drude lorentz water {float(r["parameters"]["raw_E0"]):.17g} .6 silica\n')
write('Pub_Ready/gain_dependence_linear_extended_plot.dat',[[float(x['gain_ratio']),float(x['fwhm_ps']),float(x['peak_I']),float(x['peak_time_ps']),0] for x in pulse],'NEW all-step final-cycle metrics: gain FWHM peak_I/Isat peak_time unverified_flag')
__import__('[OMEx]_fit_gain_dependence').fit_gain_dependence()
parameters='\n'.join(f'omex_Esat_{int(round(g*100))} = {e:.17g}' for g,e in scales.items())+f'\nomex_Isat_120 = {Isat:.17g}\nomex_ratio_left = {pair_config["left_eV"]/pair_config["threshold_eV"]:.17g}\nomex_ratio_right = {pair_config["right_eV"]/pair_config["threshold_eV"]:.17g}\n'
(scripts/'[OMEx]_publication_parameters.gp').write_text('# Validated fixed-Gmax scales read from NEW run metadata. No legacy scale formula.\n'+parameters)
q=scripts/'[OMEx]_paired_parameters.gp';text='# Full-integration-resolution own-peak crossing metrics; last 100000 saved-sample final mean.\n'
for j,branch in enumerate(['left','right'],1):
 x=next(x for x in pair if x['branch']==branch)
 for key,value in [('finals',x['final_last100000_saved_I']),('maxima',x['peak']),('thresholds',.81*float(x['peak'])),('t_on',x['t_on_ps'])]:text+=f'{key}[{j}] = {float(value):.17g}\n'
q.write_text(text)
last=[]
for g in [1.15,1.2,1.25]:
 a=np.loadtxt(folder(f'pulse_g{g:.2f}')/'trajectory.tsv')
 for k in range(4):
  m=(a[:,0]>=k*50)&(a[:,0]<(k+1)*50);b=a[m];write(f'Pub_Ready/sm_transition_{g:.2f}_cycle{k+1}.dat',np.c_[b[:,0]-k*50,b[:,4],b[:,5],b[:,1]*1000],'NEW time_in_cycle N Ntarget 1000*I/Isat; display sampling only')
 last.append(next(x for x in cycles if float(x['gain_ratio'])==g and float(x['cycle'])==4))
text='# New native-step pulse landmarks; intensities in 1000*I/Isat.\n'
for name,key,mult in [('peak_times','peak_time_ps',1),('peaks','peak_I',1000),('half_rises','half_rise_ps',1),('half_falls','half_fall_ps',1)]:text+=f'array {name}[3] = ['+', '.join(f'{float(x[key])*mult:.17g}' for x in last)+']\n'
(scripts/'[OMEx]_sm_transition_parameters.gp').write_text(text)
names=['pulling','probe_middle_single','freq_sweep_panels','freq_sweep_LR_panels','paired_frequency_temporal_traces','gain_switching_mechanism_1.1Gth','gain_switching_mechanism_composite','gain_dependence_linear_extended','sm_probe_peak','sm_transition']
archive=paths.IMAGES/'matplotlib_reference'
archive.mkdir(exist_ok=True)
for q in paths.IMAGES.glob('*'):
 if q.is_file() and q.suffix in ['.pdf','.png'] and not (archive/q.name).exists():shutil.copy2(q,archive/q.name)
logs=out/'logs';logs.mkdir(parents=True,exist_ok=True)
for name in names:
 q=scripts/(paths.PREFIX+name+'.gp')
 with (logs/(name+'.log')).open('w') as log:subprocess.run(['gnuplot',str(q)],cwd='/tmp',stdout=log,stderr=log,check=True)
 print('Rendered',name,flush=True)
assets=['pulling','probe_middle_single','freq_sweep_panels','freq_sweep_LR_panels','paired_frequency_temporal_traces','gain_switching_mechanism_1.1Gth','gain_switching_mechanism_composite','gain_dependence_linear_extended','sm_probe_peak_response','sm_gain_switching_transition_mechanism','sm_gain_switching_transition_cycles']
for name in assets:
 base=paths.IMAGES/(paths.PREFIX+name);subprocess.run(['pdftoppm','-png','-r','300','-singlefile',str(base)+'.pdf',str(base)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
(out/'plot_provenance.json').write_text(json.dumps({'source':'new corrected run records and full-resolution metric tables','style_source':'scripts/[OMEx]_reference/gnuplot','scripts':names,'output_directory':str(paths.IMAGES.relative_to(root)),'normalization':'validated Esat/Isat read from fixed-Gmax run metadata','legacy_seed_reused':False,'simulation_reruns':False,'png':'300 dpi raster of vector PDF','incomplete_timing':'nan, not zero','metric_tables_sha256':{n:hashlib.sha256((p/'tables'/n).read_bytes()).hexdigest() for n in ['probe_scan.csv','frequency_scan.csv','gain_dependence.csv','all_cycle_metrics.csv','paired_comparison.csv']}},indent=2))
print('PUBLICATION_STYLE_COMPLETE',len(assets),flush=True)
