from pathlib import Path
import subprocess,json
p=Path(__file__).resolve().parent;r=p.parents[2];n=p/'sources/nano_geo_matrix'
cmd=['g++','-O3','-std=c++17',*[f'-I{n/x}' for x in ['include','modules','modules/cup']],'-I/usr/include/eigen3',str(p/'driver.cxx'),'-o',str(p/'driver'),'-lgsl','-lgslcblas','-lm','-larmadillo']
(p/'build_command.json').write_text(json.dumps(cmd,indent=2));subprocess.run(cmd,check=True)
(p/'compiler.txt').write_text(subprocess.check_output(['g++','--version'],text=True))
