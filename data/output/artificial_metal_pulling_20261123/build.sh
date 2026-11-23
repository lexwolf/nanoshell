#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
library="$(realpath ../../../extern/nano_geo_matrix)"
mkdir -p overlay
cp -r "$library/modules/cup" overlay/
python3 - <<'PY'
from pathlib import Path
path = Path("overlay/cup/materials.hpp")
source = path.read_text()
anchor = "    ceps_inf = eps + ome_p2 / (ome * (ome + img * Gam_d));"
assert source.count(anchor) == 1
path.write_text(source.replace(anchor, "    eps += ::test_susceptibility(ome);\n" + anchor))
PY
git -C ../../.. rev-parse HEAD > revisions.txt
git -C "$library" rev-parse HEAD >> revisions.txt
g++ --version >> revisions.txt
mkdir -p overlay/nano_geo_matrix/quasi_static/spaser
cp "$library/include/nano_geo_matrix/quasi_static/spaser/nanoshell_intensity_steady_state.hpp" overlay/nano_geo_matrix/quasi_static/spaser/
python3 - <<'PATCH'
from pathlib import Path
path = Path("overlay/nano_geo_matrix/quasi_static/spaser/nanoshell_intensity_steady_state.hpp")
source = path.read_text()
start = source.index("double find_Omega(")
end = source.index("double gOmega(", start)
source = source[:start] + '#include "experiment_pulling.hpp"\n\n' + source[end:]
path.write_text(source)
PATCH
g++ -O3 -std=c++17 -I. -Ioverlay -I"$library/include" -I"$library/modules" -I/usr/include/eigen3 driver.cxx -o driver -lgsl -lgslcblas -larmadillo -lm
