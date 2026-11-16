"""Lossless little-endian float64 intensity storage, byte-shuffled in 65536-sample blocks.

.sf64.gz is gzip-compressed: magic b'SF64v1\\0\\0', then repeated uint32-LE
sample counts followed by 8*count bytes in byte-plane order (8,count).
No floating-point arithmetic or loss of bits occurs during packing/unpacking.
"""
import gzip,struct
from pathlib import Path
import numpy as np
MAGIC=b'SF64v1\0\0'
def iter_native(folder):
 folder=Path(folder);packed=folder/'intensity.sf64.gz'
 if packed.exists():
  with gzip.open(packed,'rb') as f:
   if f.read(8)!=MAGIC:raise ValueError('invalid shuffled intensity header')
   while h:=f.read(4):
    n=struct.unpack('<I',h)[0];raw=f.read(n*8)
    if len(raw)!=n*8:raise ValueError('truncated shuffled intensity block')
    yield np.frombuffer(raw,dtype=np.uint8).reshape(8,n).T.copy().tobytes()
 else:
  with gzip.open(folder/'intensity.f64.gz','rb') as f:
   while raw:=f.read(65536*8):yield raw
def load(folder):return np.frombuffer(b''.join(iter_native(folder)),dtype='<f8')
