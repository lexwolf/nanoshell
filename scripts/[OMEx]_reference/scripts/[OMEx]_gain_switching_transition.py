"""Diagnostic plots around 1.20 Gth; source files are never modified.
Run with python3 Gain_Switching/scripts/gain_switching_transition.py.
"""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gain_dependence as gd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'Gain_Switching/Pub_Ready'
CASES=[(1.15,ROOT/'NEW_stuff/Data/GainDependence/1.15Gth.dat'),(1.20,ROOT/'NEW_stuff/Extra/1.2Gth.dat'),(1.25,ROOT/'NEW_stuff/Extra/1.25Gth.dat')]
COLORS=['#1565C0','#D97706','#238567','#8436A8']
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':True,'axes.spines.right':True,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})

def crossing(x,y,level,up):
    hit=np.flatnonzero((y[:-1]<level)&(y[1:]>=level) if up else (y[:-1]>=level)&(y[1:]<level))
    return [float(x[i]+(level-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])) for i in hit]

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/(name+'.png'),dpi=180)
    plt.close(fig)

def main():
    fig,ax=plt.subplots(2,3,figsize=(12,6.5),sharex=True,sharey='row',layout='constrained')
    cyc,bx=plt.subplots(2,3,figsize=(12,6.5),sharex='row',sharey='row',layout='constrained')
    rows=[]; provenance=[]
    for col,(g,path) in enumerate(CASES):
        a=np.loadtxt(path); t=a[:,0]; n=a[:,6]; intensity=(a[:,5]/gd.esat(g))**2
        if not np.isfinite(a).all() or np.any(np.diff(t)<0) or t[-1]<200: raise ValueError(path)
        provenance.append(dict(gain=g,source=str(path.relative_to(ROOT)),sha256=gd.digest(path)))
        for cycle in range(4):
            mask=(t>=50*cycle)&(t<50*(cycle+1)); x=t[mask]-50*cycle; y=intensity[mask]; nn=n[mask]
            metric=gd.analyze_cycle(t,intensity,cycle*50)
            rise=crossing(x,y,y.max()/2,True); fall=crossing(x,y,y.max()/2,False)
            up=crossing(x,nn,1/g,True); down=crossing(x,nn,1/g,False)
            row=dict(gain=g,cycle=cycle+1,peak_time_ps=metric['peak_time_ps'],FWHM_ps=metric['fwhm_ps'],peak_intensity=metric['peak_intensity'],half_rise_ps=rise[0],half_fall_ps=fall[-1],N_at_peak=float(nn[np.argmax(y)]),Ncrit=1/g,N_up_ps=up[0] if up else float('nan'),N_down_ps=down[0] if down else float('nan'),N_crossing_count=len(up)+len(down))
            rows.append(row)
            # Subsampling only affects display, not numerical extraction.
            indices=np.unique(np.r_[np.arange(0,len(x),30),np.argmax(y),len(x)-1])
            bx[0,col].plot(x[indices],y[indices]*1e3,color=COLORS[cycle],lw=1.4,ls=['-','--',':','-.'][cycle],label=f'Cycle {cycle+1}')
            bx[1,col].plot(x[indices],nn[indices],color=COLORS[cycle],lw=1.4,ls=['-','--',':','-.'][cycle])
            if cycle==3:
                target=.5*(1-np.cos(2*np.pi*x/50))
                ax[0,col].plot(x[indices],target[indices],color='#D97706',ls='--',label='Prescribed inversion')
                ax[0,col].plot(x[indices],nn[indices],color='#8436A8',label='Actual inversion')
                ax[0,col].axhline(1/g,color='0.45',ls=':',label='Static threshold')
                ax[1,col].plot(x[indices],y[indices]*1e3,color='#1565C0')
                ax[1,col].fill_between(x[indices],y[indices]*1e3,color='#1565C0',alpha=.12)
                ax[1,col].hlines(y.max()*500,rise[0],fall[-1],color='#D97706',lw=2,label='FWHM')
                ax[1,col].plot(metric['peak_time_ps'],y.max()*1e3,'o',color='#1565C0',ms=4)
                ax[1,col].text(.97,.91,f"Peak: {metric['peak_time_ps']:.3f} ps\nFWHM: {metric['fwhm_ps']:.3f} ps",transform=ax[1,col].transAxes,ha='right',va='top',fontsize=9)
                for aa in ax[:,col]:
                    aa.axvline(25,color='0.6',ls='--',lw=.8)
                    for v in (up[0],down[0]): aa.axvline(v,color='#238567',ls=':',lw=1)
        for aa in (ax[0,col],bx[0,col]): aa.set_title(rf'$G_{{\max}}/G_{{\mathrm{{th}}}}={g:.2f}$')
        for aa in (ax[1,col],bx[1,col]): aa.set_xlabel('Time within cycle (ps)')
        for aa in ax[:,col]: aa.set_xlim(12,38)
        for aa in bx[:,col]: aa.set_xlim(12,38)
        ax[0,col].set_ylim(.55,1.04); ax[1,col].set_ylim(0,12.5)
        bx[0,col].set_ylim(0,12.5); bx[1,col].set_ylim(.55,1.04)
        bx[1,col].axhline(1/g,color='0.45',ls=':',lw=.8)
    ax[0,0].set_ylabel('Population inversion');ax[1,0].set_ylabel(r'$I/I_{\mathrm{sat}}$ ($\times 10^{-3}$)')
    bx[0,0].set_ylabel(r'$I/I_{\mathrm{sat}}$ ($\times 10^{-3}$)');bx[1,0].set_ylabel('Population inversion')
    ax[0,1].legend(loc='lower center',fontsize=8)
    bx[0,1].legend(loc='upper left',fontsize=8)
    fig.suptitle('Mechanism across 1.20 Gth | final cycle (150-200 ps)\nGray vertical: pump maximum; green verticals: actual-inversion threshold crossings',fontsize=11)
    cyc.suptitle('Successive cycles overlaid at equal pump phase | cycles 1-4\nCoincident curves indicate repeatability at the plotted resolution',fontsize=11)
    save(fig,'gain_switching_transition_mechanism')
    save(cyc,'gain_switching_transition_cycles')
    with (OUT/'gain_switching_transition_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'gain_switching_transition_provenance.json').write_text(json.dumps(dict(files=provenance,normalization='Esat=0.317931*sqrt(1.10/g); fixed tau1',pump='0.5*(1-cos(2*pi*t/50)), inherited from manuscript',note='Metrics use every saved sample. Display retains every 30th point plus maxima. Repeated rounded timestamps limit timing precision; no simulation or timestep convergence test performed.'),indent=2))
    for r in rows: print(r)

if __name__=='__main__': main()
