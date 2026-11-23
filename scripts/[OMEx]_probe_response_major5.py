"""Targeted continuing-probe diagnostic for Reviewer 1 major comment 5.

Actions: inspect (existing data only), refine (three half-step runs plus a dense-output control), analyze, validate.
No approved figure, manuscript, or original numerical record is overwritten.
"""
from pathlib import Path
import argparse
import ast
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import concurrent.futures
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'data/output/gain_switching_corrected_20261114'
OUT = ROOT / 'data/output/probe_response_major5_20261120'
IMG = ROOT / 'img/output/[OMEx]_probe_response_major5'
SCRIPTS = ROOT / 'scripts'
SELECTED = [('weak', 200), ('crossover', 366), ('strong', 466)]
TOL = 1e-6  # relative intensity, explicitly a diagnostic tolerance

def read_tsv(path):
    return dict(line.split('\t', 1) for line in path.read_text().splitlines())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_csv(name, rows):
    with (OUT / 'tables' / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def folder(index):
    return OLD / 'runs' / f'probe_g1.1_{index:04d}'

def native(path):
    return __import__('[OMEx]_native_io').load(path)

def inspect():
    rows = list(csv.DictReader((OLD / 'tables/probe_scan.csv').open()))
    selected = []
    evidence = []
    for regime, index in SELECTED:
        f = folder(index)
        result = json.loads((f / 'result.json').read_text())
        par = read_tsv(f / 'parameters.tsv')
        assert float(par['gain_ratio']) == 1.1 and par['pump'] == 'step'
        assert par['initial_q'] == '0,0,0' and float(par['initial_N']) == 0
        assert float(par['frequency_eV']) == float(par['threshold_frequency_eV'])
        assert par['supplied_tau1_matches'] == '1'
        seed = float(par['E0_over_Esat'])
        selected.append(dict(regime=regime, index=index, seed=seed,
                             source=str(f.relative_to(ROOT)), parameters=par,
                             configuration=result['configuration']))
    for r in rows:
        if float(r['gain_ratio']) != 1.1:
            continue
        seed = float(r['E0_over_Esat'])
        if int(r['index']) not in [200, 299, 333, 349, 366, 382, 399, 432, 466, 499]:
            continue
        times = []
        sampled_seeds = []
        for gain in [1.05, 1.1, 1.15]:
            x = min([x for x in rows if float(x['gain_ratio']) == gain],
                    key=lambda x: abs(np.log10(float(x['E0_over_Esat']) / seed)))
            assert abs(float(x['E0_over_Esat']) / seed - 1) < 1e-4
            sampled_seeds.append(float(x['E0_over_Esat']))
            times.append(float(x['t_on_ps']))
        evidence.append(dict(E0_over_Esat=seed, ton_g105_ps=times[0],
                             ton_g110_ps=times[1], ton_g115_ps=times[2],
                             gain_delay_spread_ps=max(times)-min(times),
                             actual_seed_g105=sampled_seeds[0], actual_seed_g110=sampled_seeds[1],
                             actual_seed_g115=sampled_seeds[2],
                             max_relative_seed_mismatch=max(abs(x/seed-1) for x in sampled_seeds)))
    write_csv('regime_selection.csv', evidence)
    (OUT / 'inputs/selection.json').write_text(json.dumps(selected, indent=2))
    print(json.dumps(selected, indent=2))

def refine():
    inspect()
    source = ROOT / 'src/[OMEx]_gain_switching.cxx'
    binary = OUT / 'driver.bin'
    cmd = ['g++', '-O3', '-std=c++17',
           *[f'-I{ROOT/"extern/nano_geo_matrix"/x}' for x in ['include', 'modules', 'modules/cup']],
           '-I/usr/include/eigen3', str(source), '-o', str(binary),
           '-lgsl', '-lgslcblas', '-lm', '-larmadillo']
    (OUT / 'inputs/build_command.json').write_text(json.dumps(cmd, indent=2))
    subprocess.run(cmd, check=True)
    provenance = dict(
        HEAD=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        branch=subprocess.check_output(['git', 'branch', '--show-current'], cwd=ROOT, text=True).strip(),
        library_HEAD=subprocess.check_output(['git', '-C', 'extern/nano_geo_matrix', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        initial_status=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
        compiler=subprocess.check_output(['g++', '--version'], text=True),
        source_sha256=sha(source), binary_sha256=sha(binary))
    (OUT / 'inputs/provenance.json').write_text(json.dumps(provenance, indent=2))
    shutil.copy2(source, OUT / 'inputs/driver.cxx')
    selected = json.loads((OUT / 'inputs/selection.json').read_text())
    selected.append({**selected[-1], 'kind': 'dense'})
    def run(case):
        kind = case.get('kind', 'half')
        dt = 1 if kind == 'dense' else .5
        f = OUT / 'runs' / (case['regime'] + '_' + kind)
        f.mkdir(exist_ok=True)
        config = dict(seed=case['seed'], gain=1.1, energy='threshold', pump='step',
                      T=100, dt=dt, cadence_ps=.0002, binary_sha256=sha(binary))
        result = f / 'result.json'
        if result.exists():
            previous = json.loads(result.read_text())
            assert previous['configuration'] == config and previous['returncode'] == 0
            return
        cmd = [str(binary), 'run', str(f), '1.1', 'threshold', str(case['seed']),
               'step', '100', str(dt), '.0002', '1', '0', '1', '0']
        (f / 'command.json').write_text(json.dumps(cmd, indent=2))
        with (f / 'run.log').open('w') as log:
            r = subprocess.run(cmd, stdout=log, stderr=log, timeout=120)
        result.write_text(json.dumps(dict(configuration=config, returncode=r.returncode), indent=2))
        assert r.returncode == 0
        with (f / 'intensity.f64').open('rb') as src, gzip.open(f / 'intensity.f64.gz', 'wb', compresslevel=1) as dst:
            shutil.copyfileobj(src, dst)
        (f / 'intensity.f64').unlink()
        print('Completed', case['regime'], kind, flush=True)
    with concurrent.futures.ThreadPoolExecutor(3) as ex:
        list(ex.map(run, selected))

def cross(I, dt, level):
    hits = np.flatnonzero((I[:-1] < level) & (I[1:] >= level))
    if not len(hits):
        return np.nan
    i = hits[0]
    return (i + (level-I[i])/(I[i+1]-I[i])) * dt

def stats(I, dt, lo, hi):
    a = I[int(np.ceil(lo/dt)):min(len(I), int(np.ceil(hi/dt)))]
    mean = float(np.mean(a))
    raw = float(np.ptp(a)/mean)
    floor = float(32 * np.finfo(float).eps)  # reporting floor, not a physical bound
    return dict(start_ps=lo, end_ps=min(hi, (len(I)-1)*dt), samples=len(a),
                mean_I_over_Isat=mean, min_I_over_Isat=float(a.min()), max_I_over_Isat=float(a.max()),
                measured_relative_ptp=raw, roundoff_reporting_floor=floor,
                relative_ptp_report=(f'{raw:.8g}' if raw > floor else 'unresolved_at_roundoff_floor'),
                assessment=('flat_within_tolerance' if raw < TOL else 'resolved_variation'))

def analyze():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks
    selected = json.loads((OUT / 'inputs/selection.json').read_text())
    summaries, windows, sensitivity, horizon, transients, periods = [], [], [], [], [], []
    cache = {}
    for case in selected:
        regime = case['regime']
        for kind, f in [('baseline', ROOT/case['source']), ('half', OUT/'runs'/(regime+'_half'))]:
            par = read_tsv(f/'parameters.tsv')
            dt = float(par['integration_step_ps'])
            I = native(f)
            tr = np.loadtxt(f/'trajectory.tsv')
            assert abs(float(par['E0_over_Esat'])/case['seed']-1) < 1e-12
            assert float(par['frequency_eV']) == float(par['threshold_frequency_eV'])
            assert par['pump'] == 'step' and float(par['gain_ratio']) == 1.1
            end = (len(I)-1)*dt
            peak = float(I.max()); peak_t = float(I.argmax()*dt)
            ton = cross(I, dt, .81*peak)
            rows = []
            for lo, hi in [(0,.1),(.1,.5),(.5,1),(1,2),(2,5),(5,10),(10,20),(20,40),(40,60),(60,80),(80,100)]:
                row = dict(regime=regime, integration=kind, **stats(I,dt,lo,hi))
                nt=tr[(tr[:,0]>=lo)&(tr[:,0]<hi),4]
                row['saved_trace_N_mean']=float(nt.mean())
                row['saved_trace_N_min']=float(nt.min())
                row['saved_trace_N_max']=float(nt.max())
                windows.append(row);rows.append(row)
            late = rows[-1]
            plateau = late['mean_I_over_Isat']
            overshoot = peak/plateau-1
            # No period from roundoff-level fluctuations or from Omega(th)=0.
            period_status = 'not_estimated_no_resolved_late_modulation'
            if late['assessment'] != 'flat_within_tolerance':
                period_status = 'unresolved_requires_adaptive_followup'
            periods.append(dict(regime=regime,integration=kind,status=period_status))
            early = I[:int(10/dt)]
            for prom in [1e-3, 1e-6]:
                peaks, _ = find_peaks(early, prominence=prom*plateau)
                valleys, _ = find_peaks(-early, prominence=prom*plateau)
                transients.append(dict(regime=regime,integration=kind,prominence_relative_to_late_mean=prom,
                                       local_maxima_0_10ps=len(peaks),local_minima_0_10ps=len(valleys)))
            for H in [10,20,40,60,80,100]:
                a=I[:min(len(I),int(H/dt)+1)]
                horizon.append(dict(regime=regime,integration=kind,horizon_ps=H,
                                    maximum_I=float(a.max()),own81_crossing_ps=cross(a,dt,.81*a.max())))
            last_outside = np.flatnonzero(abs(I/plateau-1) > TOL)
            settling = ((last_outside[-1]+1)*dt if len(last_outside) and last_outside[-1]+1<len(I) else (0 if not len(last_outside) else np.nan))
            late_rows = [r for r in rows if r['start_ps'] >= 40]
            mean_drift = max(abs(r['mean_I_over_Isat']/plateau-1) for r in late_rows)
            classification = ('flat_within_relative_tolerance_1e-6_over_final_windows' if all(r['assessment']=='flat_within_tolerance' for r in late_rows) and mean_drift<TOL else 'unresolved_modulation_or_slow_transient_requires_followup')
            summary = dict(regime=regime,integration=kind,E0_over_Esat=case['seed'],
                           source=str(f.relative_to(ROOT)),energy_eV=float(par['frequency_eV']),
                           raw_E0=float(par['raw_E0']),Esat=float(par['Esat']),Isat=float(par['Isat']),
                           observation_start_ps=0,observation_end_ps=end,dt_ps=dt,
                           maximum_I_over_Isat=peak,maximum_sample_time_ps=peak_t,
                           criterion_I_over_Isat=.81*peak,t_on_ps=ton,
                           late_start_ps=80,late_end_ps=end,late_mean_I_over_Isat=plateau,
                           late_measured_relative_ptp=late['measured_relative_ptp'],
                           late_relative_ptp_report=late['relative_ptp_report'],
                           late_mean_N=float(tr[tr[:,0]>=80,4].mean()),
                           peak_over_late_mean_minus_one=overshoot,
                           peak_location_class=('late_plateau_no_resolved_overshoot' if overshoot<1e-6 else 'transient_overshoot'),
                           response_class=classification,
                           finite_window_settling_time_1e_6_ps=settling,
                           max_late_window_mean_relative_drift=mean_drift,
                           late_period_status=period_status)
            summaries.append(summary)
            # Keep only modest plotting arrays; extraction above uses every integration sample.
            stride=max(1, int(.0002/dt))
            cache[(regime,kind)]=(np.arange(0,len(I),stride)*dt,I[::stride].copy(),tr,summary)
            del I
        a,b=summaries[-2:]
        sensitivity.append(dict(regime=regime,delta_ton_half_minus_baseline_ps=b['t_on_ps']-a['t_on_ps'],
                                relative_late_mean_change=b['late_mean_I_over_Isat']/a['late_mean_I_over_Isat']-1,
                                baseline_variation=a['late_relative_ptp_report'],half_variation=b['late_relative_ptp_report'],
                                baseline_class=a['response_class'],half_class=b['response_class']))
    write_csv('summary.csv',summaries);write_csv('windows.csv',windows)
    write_csv('step_sensitivity.csv',sensitivity);write_csv('horizon_checks.csv',horizon)
    write_csv('transient_extrema.csv',transients);write_csv('period_assessment.csv',periods)
    # The late zoom uses a fixed diagnostic tolerance scale, rather than magnifying roundoff into apparent oscillation.
    plt.rcParams.update({'font.family':'serif','font.size':10,'figure.constrained_layout.use':True})
    fig,axes=plt.subplots(3,3,figsize=(11.5,8.2))
    for col,(regime,_) in enumerate(SELECTED):
        t,I,tr,s=cache[(regime,'baseline')]
        th,Ih,trh,sh=cache[(regime,'half')]
        color=['#1565C0','#D97706','#8436A8'][col]
        xmax=[14,5,.65][col]
        m=t<=xmax;mh=th<=xmax
        axes[0,col].plot(t[m],I[m],color=color,label='Original step')
        axes[0,col].plot(th[mh],Ih[mh],'--',lw=.8,color='gray',label='Half step')
        axes[0,col].axhline(s['criterion_I_over_Isat'],color='gray',ls=':',label='0.81 Jmax')
        axes[0,col].plot(s['t_on_ps'],s['criterion_I_over_Isat'],'o',color=color,ms=4)
        axes[0,col].axvline(s['t_on_ps'],color=color,ls=':',lw=.8)
        axes[0,col].set(xlim=(0,xmax),ylim=(0,1.08*s['maximum_I_over_Isat']),xlabel='t (ps)',ylabel=r'$J=I/I_{sat}$',title=f'{regime.capitalize()}\n'+r'$E_0/E_{sat}=$'+f"{s['E0_over_Esat']:.6g}")
        axes[0,col].text(.97,.12,r'$t_{on}$'+f" = {s['t_on_ps']:.5g} ps",ha='right',transform=axes[0,col].transAxes,fontsize=9)
        if regime=='strong':
            tr=np.loadtxt(OUT/'runs/strong_dense/trajectory.tsv')
        mask=tr[:,0]<=xmax;maskh=trh[:,0]<=xmax
        axes[1,col].plot(tr[mask,0],tr[mask,4],color=color)
        axes[1,col].plot(trh[maskh,0],trh[maskh,4],'--',lw=.8,color='gray')
        axes[1,col].set(xlim=(0,.18 if regime=='strong' else xmax),ylim=(-.05,.25) if regime=='strong' else (0,1.03),xlabel='t (ps)',ylabel='Inversion N')
        if regime=='strong':
            axes[1,col].axhline(0,color='black',lw=.6)
            axes[1,col].text(.98,.91,'Damped inversion transient',ha='right',transform=axes[1,col].transAxes,fontsize=8.5)
        for kind,ls in [('baseline','-'),('half','--')]:
            tx,jx,_,sx=cache[(regime,kind)];mask=tx>=40
            axes[2,col].plot(tx[mask][::20],(jx[mask][::20]/sx['late_mean_I_over_Isat']-1)*100,ls,color=color,lw=.8)
        axes[2,col].set(xlim=(40,100),ylim=(-1e-4,1e-4),xlabel='t (ps)',ylabel=r'$100[J/\langle J\rangle_{80-100}-1]$ (%)')
        axes[2,col].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        axes[2,col].text(.5,.91,f"Late mean J = {s['late_mean_I_over_Isat']:.8g}\nLate mean N = {s['late_mean_N']:.7g}\nFlat within 10⁻⁶ relative tolerance",ha='center',va='top',transform=axes[2,col].transAxes,fontsize=8.5)
    axes[0,0].legend(loc='upper left',fontsize=7)
    fig.suptitle('Continuing coherent probe at exact threshold; G/Gth = 1.10, prescribed inversion = 1\nNative-step first crossings; separate intensity scales; late zooms share a stated tolerance scale',fontsize=11)
    for ext in ['pdf','png']:
        fig.savefig(IMG/f'[OMEx]_continuing_probe_responses.{ext}',dpi=200)
    plt.close(fig)
    # Save a compact, explicit validation record.
    check=dict(all_six_records_completed=True,late_relative_tolerance=TOL,
               max_abs_step_ton_difference_ps=max(abs(r['delta_ton_half_minus_baseline_ps']) for r in sensitivity),
               max_abs_step_late_mean_relative_change=max(abs(r['relative_late_mean_change']) for r in sensitivity),
               no_period_assigned_from_roundoff=True,
               note='Zero measured range in bit-identical float64 samples is labeled unresolved at the reporting floor, never physical zero modulation.')
    (OUT/'tables/validation.json').write_text(json.dumps(check,indent=2))
    dense=native(OUT/'runs/strong_dense');original=native(folder(466));bit_equal=bool(np.array_equal(dense,original));check['strong_output_cadence_native_intensity_bit_identical']=bit_equal
    check['strong_output_cadence_max_abs_intensity_difference']=float(np.max(abs(dense-original)))
    assert bit_equal, 'Output-cadence-only comparison changed native intensities'
    oldtr=np.loadtxt(folder(466)/'trajectory.tsv');newtr=np.loadtxt(OUT/'runs/strong_dense/trajectory.tsv')
    n_extrema=[]
    n_peaks=[]
    for kind in ['dense','half']:
        nt=np.loadtxt(OUT/'runs'/('strong_'+kind)/'trajectory.tsv')
        nt=nt[nt[:,0]<.5]
        up,_=find_peaks(nt[:,4],prominence=1e-4)
        down,_=find_peaks(-nt[:,4],prominence=1e-4)
        n_extrema.append(dict(integration='baseline_dense' if kind=='dense' else 'half',minimum_N=float(nt[:,4].min()),minimum_N_time_ps=float(nt[nt[:,4].argmin(),0]),maximum_N=float(nt[:,4].max()),maximum_N_time_ps=float(nt[nt[:,4].argmax(),0]),maxima_count=len(up),minima_count=len(down),N_prominence=1e-4))
        for typ,indices in [('maximum',up),('minimum',down)]:
            for i in indices:n_peaks.append(dict(integration=kind,type=typ,time_ps=float(nt[i,0]),N=float(nt[i,4])))
    write_csv('strong_inversion_extrema.csv',n_extrema)
    write_csv('strong_inversion_extrema_times.csv',n_peaks)
    write_csv('output_cadence_check.csv',[dict(case='strong',original_output_cadence_ps=float(np.median(np.diff(oldtr[:,0]))),dense_output_cadence_ps=float(np.median(np.diff(newtr[:,0]))),native_intensity_bit_identical=bit_equal,original_saved_N_max=float(oldtr[:,4].max()),dense_saved_N_max=float(newtr[:,4].max()))])
    (OUT/'tables/validation.json').write_text(json.dumps(check,indent=2))
    print(json.dumps(summaries,indent=2));print(json.dumps(sensitivity,indent=2))

def validate():
    rows=list(csv.DictReader((OUT/'tables/summary.csv').open()))
    for row in rows:
        f=ROOT/row['source']
        par=read_tsv(f/'parameters.tsv');met=read_tsv(f/'metrics.tsv')
        assert met['status']=='completed'
        assert abs(float(row['t_on_ps'])-float(met['ton_81_peak_ps']))<1e-10
        assert abs(float(par['raw_E0'])/float(par['Esat'])/float(row['E0_over_Esat'])-1)<1e-12
        assert row['late_period_status']=='not_estimated_no_resolved_late_modulation'
        assert row['response_class']=='flat_within_relative_tolerance_1e-6_over_final_windows'
    steps=list(csv.DictReader((OUT/'tables/step_sensitivity.csv').open()))
    assert max(abs(float(x['delta_ton_half_minus_baseline_ps'])) for x in steps)<1e-4
    assert max(abs(float(x['relative_late_mean_change'])) for x in steps)<1e-9
    extrema=list(csv.DictReader((OUT/'tables/strong_inversion_extrema.csv').open()))
    assert all(float(x['minimum_N'])<-.02 and int(x['maxima_count'])>=5 for x in extrema)
    assert abs(float(extrema[0]['minimum_N'])-float(extrema[1]['minimum_N']))<1e-4
    for case in ['weak','crossover','strong']:
        h=[float(x['own81_crossing_ps']) for x in csv.DictReader((OUT/'tables/horizon_checks.csv').open()) if x['regime']==case and x['integration']=='baseline' and float(x['horizon_ps'])>=20]
        assert np.ptp(h)<1e-8
    ast.parse(Path(__file__).read_text())
    print('Parameter, normalization, native-crossing, horizon, half-step and inversion-transient checks passed.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['inspect','refine','analyze','validate'])
    args=parser.parse_args()
    for d in ['inputs','runs','tables','logs']:(OUT/d).mkdir(parents=True,exist_ok=True)
    IMG.mkdir(exist_ok=True)
    globals()[args.action]()
