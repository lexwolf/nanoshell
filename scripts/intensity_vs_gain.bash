#!/bin/bash
# intensity_vs_gain.bash
# Run nanoshell_G_p3 and normalize its intensity columns by Isat.

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
range_min="0"
range_max="2"
points="201"
omega_samples="10000"
out_dir="../data/output/intensity_vs_gain"
out_file="$out_dir/intensity_vs_gain.dat"

show_help() {
    echo "Usage: bash $0 [-c] [-h] [-n points] [-r 0:max] [-k omega_samples] [-o output.dat]"
    echo ""
    echo "Options:"
    echo "  -c, --compile   Compile frohlich.cxx and nanoshell_G_p3.cxx before running"
    echo "  -h, --help      Show this help message and exit"
    echo "  -n, --points    Number of gain samples, default: 201"
    echo "  -r, --range     Gain range in G/Gth. This wrapper supports ranges starting at 0, default: 0:2"
    echo "  -k, --omega     Omega samples used inside nanoshell_G_p3, default: 10000"
    echo "  -o, --output    Output table, default: $out_file"
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
        -n|--points)
            if [ -z "${2:-}" ]; then
                echo "Error: -n|--points requires an integer argument" >&2
                exit 1
            fi
            points="$2"
            shift 2
            ;;
        -r|--range)
            if [ -z "${2:-}" ] || [[ "$2" != *:* ]]; then
                echo "Error: -r|--range requires min:max" >&2
                exit 1
            fi
            IFS=: read -r range_min range_max <<< "$2"
            shift 2
            ;;
        -k|--omega)
            if [ -z "${2:-}" ]; then
                echo "Error: -k|--omega requires an integer argument" >&2
                exit 1
            fi
            omega_samples="$2"
            shift 2
            ;;
        -o|--output)
            if [ -z "${2:-}" ]; then
                echo "Error: -o|--output requires a file path" >&2
                exit 1
            fi
            out_file="$2"
            out_dir="$(dirname "$out_file")"
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

if ! [[ "$points" =~ ^[0-9]+$ ]] || [ "$points" -lt 2 ]; then
    echo "Error: points must be an integer >= 2" >&2
    exit 1
fi

if ! [[ "$omega_samples" =~ ^[0-9]+$ ]] || [ "$omega_samples" -lt 2 ]; then
    echo "Error: omega samples must be an integer >= 2" >&2
    exit 1
fi

if [ "$range_min" != "0" ] && [ "$range_min" != "0.0" ]; then
    echo "Error: nanoshell_G_p3 sweeps from zero; use a range starting at 0, e.g. -r 0:2" >&2
    exit 1
fi

bin_targets=(../bin/fro ../bin/Gap)
src_targets=(../src/frohlich.cxx ../src/nanoshell_G_p3.cxx)

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
        echo "> Hint: rerun with the -c flag to rebuild the executables."
    fi
    exit 1
fi

if [ "$compile_requested" != true ] && [ "../src/nanoshell_G_p3.cxx" -nt "../bin/Gap" ]; then
    echo "> Source is newer than ../bin/Gap; recompiling."
    compile_requested=true
fi

if [ "$compile_requested" = true ]; then
    echo
    echo "> Compiling codes..."
    $CXX $CXXFLAGS ../src/frohlich.cxx -o ../bin/fro $LDFLAGS $LIBS
    $CXX $CXXFLAGS ../src/nanoshell_G_p3.cxx -o ../bin/Gap $LDFLAGS $LIBS
    echo "> ...Done!"
    echo
fi

file_path="../data/input/nanosphere_eV.dat"
if [ ! -f "$file_path" ]; then
    echo "Error: input file not found: $file_path" >&2
    exit 1
fi

fro_output=( $("../bin/fro") )
omega_th="${fro_output[0]:-}"
Gth="${fro_output[1]:-}"

if [ -z "$omega_th" ] || [ -z "$Gth" ]; then
    echo "Error: frohlich calculation failed (omega_th or Gth empty)." >&2
    exit 1
fi

gain_steps="$((points - 1))"
raw_file="../data/output/emission_maximum.dat"
warn_file="$out_dir/intensity_vs_gain.warnings.log"

mkdir -p "$out_dir"
: > "$warn_file"

echo "> Frohlich frequency omega_th = $omega_th eV"
echo "> Threshold gain Gth = $Gth"
echo "> Running nanoshell_G_p3 for G/Gth = 0 ... $range_max ($points points)"
echo

"../bin/Gap" "$range_max" "$gain_steps" "$omega_samples" > "$warn_file" 2>&1

if [ ! -f "$raw_file" ]; then
    echo "Error: expected raw output not found: $raw_file" >&2
    exit 1
fi

{
    echo "# Normalized dipole emission intensity vs gain"
    echo "# generated from ../data/output/emission_maximum.dat by $0"
    echo "# omega_th = $omega_th"
    echo "# Gth = $Gth"
    echo "# raw columns: G G_over_Gth abs_p_ss_sq abs_p_analytic_sq abs_p_num_sq emission_width kex1 Isat"
    echo "# columns: G G_over_Gth abs_p_ss_sq_over_Isat abs_p_num_sq_over_Isat abs_p_analytic_sq_over_Isat Isat"
    awk '
        NF >= 8 && $1 !~ /^#/ {
            singular = ($8 == 0 || $8 == "inf" || $8 == "Inf")
            pss = (singular ? 0 : $3 / $8)
            pnum = (singular ? 0 : $5 / $8)
            panl = ($2 >= 1 && !singular ? $4 / $8 : "NaN")
            printf "%.12g %.12g %.12g %.12g %s %s\n", $1, $2, pss, pnum, panl, $8
        }
    ' "$raw_file"
} > "$out_file"

echo "> Wrote normalized intensity comparison table to $out_file"
echo "> Raw nanoshell_G_p3 table: $raw_file"
if [ -s "$warn_file" ]; then
    echo "> Wrote nanoshell_G_p3 output/warnings to $warn_file"
fi
echo "> Done."
