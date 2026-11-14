import batch as b,json,concurrent.futures,time
p=b.p; b.jobs=[];G=b.G
for g,gr in G.items():
 old={}
 for row in gr['frequency_source_rows']:old[row[1]]=max(old.get(row[1],0),row[0])
 for i,e in enumerate(gr['frequencies']):
  t=old[e];T=min(2000,max(400,1.5*t+220));minimum=1.2*t+40
  if abs(e-2.8122)<1e-10:T=6300;minimum=40
  for factor in [.5]:b.add(f'residual_tight_g{g}_{i:04d}',g,e,1e-9,'step',T,factor=factor,diag=1,minobs=minimum)
  # Match observation budgets for incomplete baseline cases; retain complete default runs.
  folder=p/'runs'/f'freq_g{g}_{i:04d}';r=folder/'result.json'
  completed=r.exists() and json.loads(r.read_text())['metrics'].get('status')=='completed'
  if not completed or abs(e-2.8122)<1e-10:b.add(f'residual_default_g{g}_{i:04d}',g,e,1e-9,'step',T,diag=1,minobs=minimum)
(p/'inputs/residual_jobs.json').write_text(json.dumps(b.jobs,indent=2))
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(16) as ex:res=list(ex.map(b.run,b.jobs))
 (p/'tables/residual_rerun_results.json').write_text(json.dumps(res,indent=2));print('RESIDUAL_BATCH_FINISHED',len(res),flush=True)
