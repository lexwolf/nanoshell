#!/bin/bash
# Supplementary Fig. 1 of [nanoph-2024-0491]: seed-field sweep.
# Uses the archived raw seed fields; published legend ratios are approximate.
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

show_help() {
    echo "Usage: bash $0 [-c|--compile] [---paper-values] [--plot|--plot-only] [-h|--help]"
    echo "Generate four time traces at 2.811 eV and G/Gth = 1.01."
    echo "Raw E0 = 1e-8, 1e-9, 1e-10, 1e-11, as in the paper archive."
    echo "Inputs: data/input/nanosphere_eV.dat and data/input/time.dat."
    echo "Data: data/output/[nanoph-2024-0491]_rangeE0/."
    echo "  -c, --compile  Build fro, Esat, and nsh before running."
    echo "  ---paper-values  Read historical inputs from $HOME/Nextcloud/Work/Research/Papers/Archive/ns_paper/versions/Draft 2/img/work/Supp/in."
    echo "                   Use archived T and pump time directly (already in ps)."
    echo "                   Both repo input files are restored on exit."
    echo "  --plot         Render the PDF after completing the sweep."
    echo "  --plot-only    Render an existing completed sweep without simulation."
    echo "  -h, --help     Show this help message."
}

compile_requested=false
paper_values=false
plot_requested=false
plot_only=false
work_dir="../data/output/[nanoph-2024-0491]_rangeE0"
plot_file="[nanoph-2024-0491]_rangeE0.gp"
render_plot() {
    mkdir -p ../img/output
    gnuplot "$plot_file"
    echo "> Rendered ../img/output/[nanoph-2024-0491]_rangeE0.pdf"
}
while [ "$#" -gt 0 ]; do
    case "$1" in
        -c|--compile) compile_requested=true ;;
        ---paper-values) paper_values=true ;;
        --plot) plot_requested=true ;;
        --plot-only) plot_only=true ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Error: unexpected argument: $1" >&2; exit 1 ;;
    esac
    shift
done
if "$plot_requested" || "$plot_only"; then
    command -v gnuplot >/dev/null || { echo "Error: gnuplot is required." >&2; exit 1; }
fi
if "$plot_only"; then
    render_plot
    exit 0
fi

if "$compile_requested"; then
    ngm_root="$(cd ../extern/nano_geo_matrix && pwd)"
    includes=(-Wall "-I$ngm_root/include" "-I$ngm_root/modules"
              "-I$ngm_root/modules/cup" -I/usr/local/include -I/usr/include/eigen3)
    sources=(frohlich Esat nanoshell_num_anl)
    targets=(fro Esat nsh)
    mkdir -p ../bin
    echo "> Compiling codes..."
    for i in "${!sources[@]}"; do
        "${CXX:-g++}" "${includes[@]}" "../src/${sources[i]}.cxx" \
            -o "../bin/${targets[i]}" -L/usr/local/lib -lgsl -lgslcblas -lm -larmadillo
    done
fi

for executable in fro Esat nsh; do
    if [ ! -x "../bin/$executable" ]; then
        echo "Error: missing executable ../bin/$executable; rerun with --compile." >&2
        exit 1
    fi
done

file_path="../data/input/nanosphere_eV.dat"
time_path="../data/input/time.dat"
for input in "$file_path" "$time_path"; do
    [ -r "$input" ] || { echo "Error: cannot read $input" >&2; exit 1; }
done
paper_input_dir="$HOME/Nextcloud/Work/Research/Papers/Archive/ns_paper/versions/Draft 2/img/work/Supp/in"
if "$paper_values"; then
    for name in nanosphere_eV.dat time.dat; do
        [ -r "$paper_input_dir/$name" ] || {
            echo "Error: cannot read historical input $paper_input_dir/$name" >&2
            exit 1
        }
    done
fi

gain_factor=1.01
omega=2.811
mkdir -p "$work_dir"
# Prevent plotting stale metadata if this run fails partway through.
rm -f "$work_dir/parameters.gp"

# The executables use fixed repo input paths. Restore the original bytes on exit.
backup_dir="$(mktemp -d)"
cp "$file_path" "$backup_dir/nanosphere_eV.dat"
cp "$time_path" "$backup_dir/time.dat"
cleanup() {
    cp "$backup_dir/nanosphere_eV.dat" "$file_path"
    cp "$backup_dir/time.dat" "$time_path"
    rm -rf "$backup_dir"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

if "$paper_values"; then
    read -r historical_T historical_pump \
        < <(awk 'NF && $1 !~ /^#/ {print; exit}' "$paper_input_dir/time.dat")
    if ! awk -v t="${historical_T:-}" -v p="${historical_pump:-}" \
        'BEGIN {exit !(t+0 > 0 && p != "" && p+0 >= 0 && p+0 <= t+0)}'; then
        echo "Error: invalid historical simulation or pump time." >&2
        exit 1
    fi
    cp "$paper_input_dir/nanosphere_eV.dat" "$file_path"
    # This paper archive already uses Ep/hbar, just like the current solver.
    cp "$paper_input_dir/time.dat" "$time_path"
    echo "> Using paper archive inputs from $paper_input_dir"
    echo "> Simulation time (ps): T = $historical_T; pump = $historical_pump"
fi

read -r a Dome ome_g G omemi omema metal model active solvent E0 rho host \
    < <(awk 'NF && $1 !~ /^#/ {print; exit}' "$file_path")
[ -n "${host:-}" ] || { echo "Error: expected 13 input parameters." >&2; exit 1; }
cp "$file_path" "$work_dir/input.dat"
cp "$time_path" "$work_dir/time.dat"

fro_result="$(../bin/fro)"
read -r omega_th Gth <<< "$fro_result"
if ! awk -v w="${omega_th:-}" -v g="${Gth:-}" \
    'BEGIN {exit !(w+0 > 0 && g+0 != 0)}'; then
    echo "Error: invalid threshold output: $fro_result" >&2
    exit 1
fi
new_G="$(awk -v f="$gain_factor" -v g="$Gth" 'BEGIN {printf "%.15g", f*g}')"

# Esat calculates its own threshold and accepts G/Gth as its argument.
esat_result="$(../bin/Esat "$gain_factor")"
esat="$(awk '$2 == "nEsat" && $3 == "=" {print $4}' <<< "$esat_result")"
if ! awk -v e="$esat" 'BEGIN {exit !(e+0 > 0)}'; then
    echo "Error: invalid saturation-field output: $esat_result" >&2
    exit 1
fi
Isat="$(awk -v e="$esat" 'BEGIN {printf "%.15g", e*e}')"
if "$paper_values"; then
    # Exact intensity normalization from the archived dull_rangeE0.gp.
    Isat=0.110087
fi
read -r T tpump < <(awk 'NF && $1 !~ /^#/ {print; exit}' "$time_path")
{
    echo "# Run parameters; saturation convention: src/Esat.cxx."
    echo "# Input snapshots: input.dat and time.dat in this directory."
    echo "paper_values = $($paper_values && echo 1 || echo 0)"
    echo "gain_factor = $gain_factor"
    echo "Gth = $Gth"
    echo "G = $new_G"
    echo "omega_th = $omega_th"
    echo "omega = $omega"
    echo "Esat = $esat"
    echo "Isat = $Isat"
    echo "raw_seed_sweep = 1"
    echo "T = $T"
    echo "tpump = $tpump"
} > "$work_dir/parameters.gp.pending"

echo "> omega_th = $omega_th eV; Gth = $Gth; G = $new_G"
echo "> Esat = $esat (src/Esat.cxx convention)"
echo "> Isat = $Isat"
echo "> Calculating for the archived raw seed fields..."
for exponent in 8 9 10 11; do
    new_E0="1.e-$exponent"
    {
        echo "$a $Dome $omega_th $new_G $omemi $omema $metal $model $active $solvent $new_E0 $rho $host"
        echo "# a Dome omega_th G omemi omema metal model active solvent E0 rho host"
    } > "$file_path"
    ratio="$(awk -v e="$new_E0" -v sat="$esat" 'BEGIN {printf "%.8g", e/sat}')"
    echo "> E0 = $new_E0; E0/Esat = $ratio (published legend is approximate)"
    ../bin/nsh "$omega"
    cp ../data/output/numtime.dat "$work_dir/numtime$(printf '%02d' "$exponent").dat"
    cp "$file_path" "$work_dir/input$(printf '%02d' "$exponent").dat"
done
mv "$work_dir/parameters.gp.pending" "$work_dir/parameters.gp"
if "$plot_requested"; then
    render_plot
fi
echo "> Done. Data saved to $work_dir; both original inputs restored on exit."
