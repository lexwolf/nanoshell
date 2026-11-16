from pathlib import Path
import json,csv,numpy as np
paths=__import__("[OMEx]_paths");p=paths.DATA
results=json.loads((p/'tables/all_results.json').read_text())+json.loads((p/'tables/residual_rerun_results.json').read_text());errors=[]
for r in results:
 a=r['parameters'];c=r['configuration'];
 if r['returncode']!=0:errors.append((c['name'],'returncode',r['returncode']));continue
 if not np.isclose(float(a['raw_E0'])/float(a['Esat']),c['seed'],rtol=1e-12,atol=0):errors.append((c['name'],'seed mismatch'))
 if not np.isclose(float(a['Isat']),float(a['Esat'])**2,rtol=1e-12):errors.append((c['name'],'square mismatch'))
 if a['supplied_tau1_matches']!='1' or a['tau1_policy']!='require_match':errors.append((c['name'],'relaxation validation'))
 if float(a['Gmax'])<=0 or a['gain_convention']!='Gmax=fixed_concentration_full_inversion':errors.append((c['name'],'gain convention'))
# Recompute pulse metrics independently from all-step native records of representative cases.
load=__import__("[OMEx]_native_io").load
for g in [1.05,1.1,1.15,1.2,1.25,1.5]:
 folder=p/f'runs/pulse_g{g:.2f}';a=dict(l.split('\t',1) for l in (folder/'parameters.tsv').read_text().splitlines());step=float(a['integration_step_ps']);I=load(folder);start=int(np.ceil(150/step));end=min(len(I),int(np.ceil(200/step)));x=I[start:end];k=int(np.argmax(x));peak=x[k];half=peak/2;ris=np.flatnonzero((x[:-1]<half)&(x[1:]>=half));fal=np.flatnonzero((x[:-1]>=half)&(x[1:]<half));ris=ris[ris<k][-1];fal=fal[fal>=k][0]
 def cross(j):return ((start+j)+(half-x[j])/(x[j+1]-x[j]))*step-150
 cyc=np.genfromtxt(folder/'cycles.tsv',names=True,delimiter='\t')[-1]
 if not np.isclose(cross(fal)-cross(ris),cyc['fwhm_ps'],rtol=1e-11,atol=1e-11):errors.append((str(g),'independent width mismatch'))
 if peak!=cyc['peak_I']:errors.append((str(g),'independent peak mismatch'))
# Periodic target must be the explicit nonzero-frequency waveform.
d=np.loadtxt(p/'runs/pulse_g1.10/trajectory.tsv');pump_error=float(np.max(np.abs(d[:,5]-.5*(1-np.cos(2*np.pi*d[:,0]/50)))))
if pump_error>1e-12:errors.append(('pump',pump_error))
a=np.fromfile(p/'pilot/strong/intensity.f64');b=np.fromfile(p/'pilot/strong_half/intensity.f64');strong=float(np.max(np.abs(a-b[::2][:len(a)]))/a.max())
out={'runs_audited':len(results),'normalization_mapping_errors':len(errors),'errors':errors,'periodic_pump_max_abs_error':pump_error,'strong_probe_dt_halving_max_difference_relative_to_peak':strong,'native_metric_independent_reconstruction':'passed' if not errors else 'failed','full_resolution':'float64 each integration step; lossless storage','time_reference':'legacy pre-step timestamp retained for post-step state'}
(p/'tables/validation.json').write_text(json.dumps(out,indent=2));print(out)
if errors:raise SystemExit(1)
