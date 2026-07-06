#!/bin/bash
# [nanoph-2024-0491]_Grange.bash
# Reconstructs data and renders the Grange figure used in [nanoph-2024-0491].

set -euo pipefail
export LC_NUMERIC="en_US.UTF-8"

script_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$script_dir" || exit 1

CXX=g++
NGM_ROOT="$(realpath ../extern/nano_geo_matrix)"
CXXFLAGS="-Wall -I$NGM_ROOT/include -I$NGM_ROOT/modules -I$NGM_ROOT/modules/cup -I/usr/local/include -I/usr/include/eigen3"
LDFLAGS="-L/usr/local/lib"
LIBS="-lgsl -lgslcblas -lm -larmadillo"

compile_requested=false
plot_requested=false
work_dir="../data/output/Grange"
gnuplot_file="[nanoph-2024-0491]_Grange.gp"
gain_factors=(1.01 1.25 1.5 1.75 2)
gain_max="2"
gap_steps="500"
gap_omega_samples="10000"

grange_a="10"
grange_Dome=".15"
grange_ome_g="2.8121972"
grange_omemi="2."
grange_omema="3.4"
grange_metal="silver"
grange_model="drude"
grange_gain_model="lorentz"
grange_solvent="water"
grange_E0="1.e-8"
grange_rho="0.6"
grange_host="silica"
grange_time_T="600"
grange_time_pump="10"

show_help() {
    echo "Usage: bash $0 [-c] [-h] [-w work_dir] [--plot]"
    echo ""
    echo "Options:"
    echo "  -c, --compile       Compile required executables before running"
    echo "  -h, --help          Show this help message and exit"
    echo "  -w, --work-dir      Output data root, default: $work_dir"
    echo "  --plot              Render the PDF after generating the data"
    echo "  --gains list        Space-separated G/Gth values, default: ${gain_factors[*]}"
    echo "  --gap-steps n       Gain steps for bin/Gap, default: $gap_steps"
    echo "  --gap-omega n       Omega samples for bin/Gap, default: $gap_omega_samples"
}

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--compile)
            compile_requested=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        -w|--work-dir)
            work_dir="$2"
            shift 2
            ;;
        --plot)
            plot_requested=true
            shift
            ;;
        --gains)
            IFS=' ' read -r -a gain_factors <<< "$2"
            shift 2
            ;;
        --gap-steps)
            gap_steps="$2"
            shift 2
            ;;
        --gap-omega)
            gap_omega_samples="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            show_help
            exit 1
            ;;
        *)
            echo "Unexpected argument: $1" >&2
            show_help
            exit 1
            ;;
    esac
done

if [ "$plot_requested" = true ] && [ ! -f "$gnuplot_file" ]; then
    echo "Error: gnuplot file not found: $gnuplot_file" >&2
    exit 1
fi

bin_targets=(../bin/fro ../bin/Gap ../bin/oap)
src_targets=(../src/frohlich.cxx ../src/nanoshell_G_p3.cxx ../src/nanoshell_ome_al_p3.cxx)

missing=()
missing_exec=false

if [ "$compile_requested" = true ]; then
    for src in "${src_targets[@]}"; do
        [ -f "$src" ] || missing+=("$src (source)")
    done
else
    for bin in "${bin_targets[@]}"; do
        if [ ! -x "$bin" ]; then
            missing_exec=true
            missing+=("$bin (executable)")
        fi
    done
fi

if [ ${#missing[@]} -gt 0 ]; then
    if [ "$compile_requested" = true ]; then
        echo "> Compilation requested but the following files are missing:"
    else
        echo "> Missing required files:"
    fi
    for item in "${missing[@]}"; do
        echo "  - $item"
    done
    if [ "$compile_requested" != true ] && [ "$missing_exec" = true ]; then
        echo "> Hint: rerun with -c to rebuild the executables."
    fi
    exit 1
fi

if [ "$compile_requested" = true ]; then
    echo "> Compiling required executables..."
    $CXX $CXXFLAGS ../src/frohlich.cxx -o ../bin/fro $LDFLAGS $LIBS
    $CXX $CXXFLAGS ../src/nanoshell_G_p3.cxx -o ../bin/Gap $LDFLAGS $LIBS
    $CXX $CXXFLAGS ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap $LDFLAGS $LIBS
fi

repo_input="../data/input/nanosphere_eV.dat"
repo_time="../data/input/time.dat"
repo_input_bak="$(mktemp)"
repo_time_bak="$(mktemp)"
cp -f "$repo_input" "$repo_input_bak"
if [ -f "$repo_time" ]; then
    cp -f "$repo_time" "$repo_time_bak"
else
    : > "$repo_time_bak"
fi
cleanup() {
    mv -f "$repo_input_bak" "$repo_input"
    if [ -s "$repo_time_bak" ]; then
        mv -f "$repo_time_bak" "$repo_time"
    else
        rm -f "$repo_time"
        rm -f "$repo_time_bak"
    fi
}
trap cleanup EXIT INT TERM

mkdir -p "$work_dir/in" "$work_dir/logs" "../data/output/oGp" "../img/output"

write_grange_input() {
    local omega="$1"
    local G_val="$2"
    {
        echo "$grange_a $grange_Dome $omega $G_val $grange_omemi $grange_omema $grange_metal $grange_model $grange_gain_model $grange_solvent $grange_E0 $grange_rho $grange_host"
        echo "# a Dome omega_th G omemi omema metal model gain_model solvent E0 rho host"
    } > "$repo_input"
}

write_grange_time() {
    {
        echo "$grange_time_T          $grange_time_pump"
        echo "#T          t_pump"
    } > "$repo_time"
}

write_grange_input "$grange_ome_g" "0."
write_grange_time
cp -f "$repo_input" "$work_dir/in/nanosphere_eV.dat"
cp -f "$repo_time" "$work_dir/in/time.dat"

read -r a Dome ome_g G omemi omema metal model gain_model solvent E0 rho host < "$repo_input"

fro_output=( $("../bin/fro") )
omega_th="${fro_output[0]:-}"
Gth="${fro_output[1]:-}"
if [ -z "$omega_th" ] || [ -z "$Gth" ]; then
    echo "Error: frohlich calculation failed." >&2
    exit 1
fi

write_input() {
    local G_val="$1"
    write_grange_input "$omega_th" "$G_val"
}

tag_for_factor() {
    local factor="$1"
    printf "%s" "$factor" | sed 's/\.//g'
}

folder_for_factor() {
    local factor="$1"
    printf "oG%sp" "$factor"
}

echo "> omega_th = $omega_th"
echo "> Gth      = $Gth"

echo "> Running global gain sweep with bin/Gap..."
"../bin/Gap" "$gain_max" "$gap_steps" "$gap_omega_samples" > "$work_dir/logs/Gap.log" 2>&1
cp -f "../data/output/emission_maximum.dat" "$work_dir/emission_maximum.dat"
cp -f "../data/output/emission_maximum.dat" "$work_dir/emission_maximum.0.dat"

for factor in "${gain_factors[@]}"; do
    folder="$(folder_for_factor "$factor")"
    factor_tag="$(tag_for_factor "$factor")"
    fixed_dir="$work_dir/$folder"
    mkdir -p "$fixed_dir"

    G_val="$(awk -v f="$factor" -v gth="$Gth" 'BEGIN { printf "%.12g", f*gth }')"
    echo "> Running fixed spectrum for G = $factor Gth = $G_val"
    write_input "$G_val"
    cp -f "$repo_input" "$work_dir/in/nanosphere_eVG_${factor}.dat"
    cp -f "$repo_input" "$work_dir/in/nanosphere_eVG${factor}.dat"
    cp -f "$repo_input" "$work_dir/in/nanosphere_eVG_${factor_tag}.dat"

    "../bin/oap" > "$fixed_dir/oap.log" 2>&1

    cp -f "../data/output/oGp/ome_p3.dat" "$fixed_dir/ome_p3.dat"
    cp -f "../data/output/oGp/ome_al.dat" "$fixed_dir/ome_al.dat"
    cp -f "../data/output/oGp/ome_ex.dat" "$fixed_dir/ome_ex.dat"
done

if [ "$plot_requested" = true ]; then
    echo "> Rendering ../img/output/[nanoph-2024-0491]_Grange.pdf"
    gnuplot "$gnuplot_file"
fi

echo "> Wrote Grange data tree to $work_dir"
echo "> Done. Original repo input restored automatically."
