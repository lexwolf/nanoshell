"""Descriptive fits to the existing corrected pulse peaks; no dynamics are rerun."""
import csv,json
import numpy as np
paths=__import__('[OMEx]_paths')
def fit_gain_dependence():
 rows=list(csv.DictReader((paths.DATA/'tables/gain_dependence.csv').open()))
 g=np.array([float(r['gain_ratio']) for r in rows]);I=np.array([float(r['peak_I']) for r in rows]);order=np.argsort(g);g=g[order];I=I[order]
 low=(g>=1.05-1e-12)&(g<=1.09+1e-12);high=g>=1.25-1e-12
 x=g[low]-1.05;y=np.log10(I[low]);m,c=np.polyfit(x,y,1);res=y-(m*x+c);R2=1-float(res@res)/float(((y-y.mean())**2).sum());a,b=np.polyfit(g[high],I[high],1);r=I[high]-(a*g[high]+b)
 local=np.diff(y)/np.diff(g[low]);report={'low_gain':{'gain_window':[1.05,1.09],'points':int(low.sum()),'fit':'log10(Ipeak/Isat)=m*(g-1.05)+c; g=Gmax/Gth','m_decades_per_unit_gain_ratio':float(m),'c':float(c),'alpha_natural_exponential_per_unit_gain_ratio':float(m*np.log(10)),'factor_per_gain_increment_0p01':float(10**(m*.01)),'R_squared_log10':R2,'max_abs_log10_residual':float(abs(res).max()),'local_segment_slopes_decades_per_unit_gain_ratio':local.tolist(),'qualification':'Descriptive approximation over five deterministic points. Curvature remains; not a demonstrated exact exponential law or statistical uncertainty estimate.'},'high_gain':{'gain_window':[1.25,float(g.max())],'points':int(high.sum()),'fit':'Ipeak/Isat=a*g+b','a':float(a),'b':float(b),'R_squared_linear':1-float(r@r)/float(((I[high]-I[high].mean())**2).sum()),'qualification':'Descriptive fit on the specified high-gain interval; not an integration-error bound.'},'interpretation':'Rapid approximately exponential rise at low gain, followed by a crossover toward the nearly linear high-gain dependence; no second exponential regime is inferred.'}
 (paths.DATA/'tables/gain_dependence_fits.json').write_text(json.dumps(report,indent=2))
 with (paths.DATA/'tables/low_gain_fit_points.csv').open('w',newline='') as f:
  w=csv.writer(f,lineterminator="\n");w.writerow(['gain_ratio','Ipeak_over_Isat','log10_I','fitted_log10_I','log10_residual']);w.writerows(zip(g[low],I[low],y,m*x+c,res))
 text=f'''# Descriptive fits generated from corrected native-step pulse peak summaries.
# m is d(log10(Ipeak/Isat))/d(Gmax/Gth), fitted over 1.05--1.09.
omex_low_m = {m:.17g}
omex_low_c = {c:.17g}
omex_low_alpha = {m*np.log(10):.17g}
omex_high_a = {a:.17g}
omex_high_b = {b:.17g}
omex_low_fit(x) = (x>=1.05 ? 10**(omex_low_m*(x-1.05)+omex_low_c) : 1/0)
'''
 (paths.ROOT/'scripts/[OMEx]_gain_dependence_fits.gp').write_text(text)
 return report
if __name__=='__main__':print(json.dumps(fit_gain_dependence(),indent=2))
