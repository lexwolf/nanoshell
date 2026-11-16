from pathlib import Path
import subprocess,json
paths=__import__("[OMEx]_paths");p=paths.DATA;r=p.parents[2];n=p/'sources/nano_geo_matrix'
paths.BINARY.parent.mkdir(parents=True,exist_ok=True)
cmd=['g++','-O3','-std=c++17',*[f'-I{n/x}' for x in ['include','modules','modules/cup']],'-I/usr/include/eigen3',str(paths.DRIVER_SOURCE),'-o',str(paths.BINARY),'-lgsl','-lgslcblas','-lm','-larmadillo']
(p/'build_command_current.json').write_text(json.dumps(cmd,indent=2));subprocess.run(cmd,check=True)
(p/'compiler_current.txt').write_text(subprocess.check_output(['g++','--version'],text=True))
