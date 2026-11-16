"""Update full-resolution metrics/diagnostics, then draw publication assets with gnuplot."""
from pathlib import Path
import subprocess,sys
scripts=Path(__file__).resolve().parent
for name in ['extract_metrics','publication_plots']:
 subprocess.run([sys.executable,str(scripts/('[OMEx]_'+name+'.py'))],check=True)
