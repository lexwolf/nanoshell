#!/bin/bash
# Validate the shared detector against saved paper inputs; preserve all paper files.
set -euo pipefail
export LC_ALL=C
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
ngm="$(cd ../extern/nano_geo_matrix && pwd)"
source_dir='../data/output/[nanoph-2024-0491]_rangeE0'
output_dir='../data/output/convergence_validation'
if [ -e "$output_dir" ]; then
    echo "Refusing to overwrite existing validation directory: $output_dir" >&2
    exit 1
fi
for n in 08 09 10 11; do
    test -r "$source_dir/input$n.dat"
    test -r "$source_dir/numtime$n.dat"
done
test -r "$source_dir/time.dat"
mkdir -p ../bin
"${CXX:-g++}" -std=c++17 -O2 -Wall -Wextra -I"$ngm/include" -I"$ngm/modules" -I"$ngm/modules/cup" \
    -I/usr/include/eigen3 ../src/emission_convergence.cxx -o ../bin/emission_convergence \
    -lgsl -lgslcblas -lm -larmadillo
mkdir -p "$output_dir/inputs"
cp "$source_dir"/input{08,09,10,11}.dat "$source_dir/time.dat" "$source_dir/parameters.gp" "$output_dir/inputs/"
git -C "$ngm" rev-parse HEAD > "$output_dir/dependency_commit.txt"
for mode in fixed automatic tight; do
    factor=1
    solver_mode="$mode"
    if [ "$mode" = tight ]; then factor=.5; solver_mode=automatic; fi
    for n in 08 09 10 11; do
        out="$output_dir/$mode/$n"
        mkdir -p "$out"
        echo "> $mode: seed 1e-$((10#$n))"
        ../bin/emission_convergence "$output_dir/inputs/input$n.dat" "$output_dir/inputs/time.dat" \
            "$out" "$solver_mode" "$factor" > "$out/run.log" 2>&1
    done
done
python3 "summarize_convergence.py"
mkdir -p ../img/output
gnuplot "convergence_validation.gp"
echo "> Completed validation in $output_dir"
