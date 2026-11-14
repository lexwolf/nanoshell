from pathlib import Path
import gzip,struct,hashlib,json,concurrent.futures,time
import numpy as np
from native_io import MAGIC,iter_native
p=Path(__file__).resolve().parent

def pack(q):
 folder=q.parent;target=folder/'intensity.sf64.gz';oldbytes=q.stat().st_size
 if target.exists():return {'folder':folder.name,'already_packed':True}
 digest=hashlib.sha256();count=0
 with gzip.open(q,'rb') as src,gzip.open(target,'wb',compresslevel=4) as dst:
  dst.write(MAGIC)
  while raw:=src.read(65536*8):
   if len(raw)%8:raise ValueError('misaligned native intensity file')
   digest.update(raw);n=len(raw)//8;count+=n;dst.write(struct.pack('<I',n));dst.write(np.frombuffer(raw,dtype=np.uint8).reshape(n,8).T.copy().tobytes())
 check=hashlib.sha256()
 for raw in iter_native(folder):check.update(raw)
 if check.digest()!=digest.digest():raise ValueError('intensity lossless roundtrip failed')
 meta={'folder':folder.name,'samples':count,'uncompressed_sha256':digest.hexdigest(),'gzip_before_bytes':oldbytes,'shuffled_gzip_bytes':target.stat().st_size,'roundtrip_verified':True};(folder/'intensity_storage.json').write_text(json.dumps(meta,indent=2));q.unlink();return meta
if __name__=='__main__':
 files=list((p/'runs').glob('*/intensity.f64.gz'));t=time.monotonic()
 existing=[json.loads(q.read_text()) for q in (p/'runs').glob('*/intensity_storage.json')]
 with concurrent.futures.ThreadPoolExecutor(8) as ex:
  rows=[]
  for i,r in enumerate(ex.map(pack,files)):
   rows.append(r)
   if (i+1)%300==0:print(i+1,'packed',time.monotonic()-t,flush=True)
 byname={r['folder']:r for r in existing+rows};(p/'tables/storage_manifest.json').write_text(json.dumps(list(byname.values()),indent=2));print('PACKING_DONE',len(rows),time.monotonic()-t,flush=True)
