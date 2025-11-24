#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"

script_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$script_dir" || exit 1

compile_requested=false
range_override=""
rho_range_override=""

show_help() {
    echo "Usage: bash $0 [-c] [-h] [-r omi:oma] [-rho \"0.4 0.5 0.6 0.7 0.8\"]"
    echo ""
    echo "Options:"
    echo "  -c    --compile    Compile the codes before executing the script"
    echo "  -h    --help       Show this help message and exit"
    echo "  -r    --range      Override ome_min:ome_max (eV) e.g. 2.0:3.5"
    echo "  -rho  --rho-range  Space-separated, increasing positive numbers for rho values (default: 0.4 0.5 0.6 0.7 0.8)"
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
        -r|--range)
            if [ -z "${2:-}" ]; then
                echo "Error: -r|--range requires an argument in the form omemi:omema" >&2
                exit 1
            fi
            range_override="$2"
            shift 2
            ;;
        -rho|--rho-range)
            if [ -z "${2:-}" ]; then
                echo "Error: -rho|--rho-range requires a space-separated list of positive numbers" >&2
                exit 1
            fi
            rho_range_override="$2"
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

if [ -n "$range_override" ] && [[ "$range_override" != *:* ]]; then
    echo "Error: -r|--range requires omemi:omema" >&2
    exit 1
fi

bin_targets=(../bin/rho2ome_sp ../bin/oap)
src_targets=(../src/rho2ome_sp.cxx ../src/nanoshell_ome_al_p3.cxx)
gp_scripts=("[nanoph-2024-0491]_rho_emi.gp" "nano-shell-sketch.bash")
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

for gp in "${gp_scripts[@]}"; do
    [ -f "$gp" ] || missing+=("$gp (script)")
done

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

if [ "$compile_requested" = true ]; then
    echo
    echo "> Compiling codes..."
    g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
    g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo
    echo "> ...Done!"
    echo
fi

none=1.e-30
file_path="../data/input/nanosphere_eV.dat"
echo "> WARNING: $file_path will be temporarily overwritten during the calculations."

read  a  Dome  ome21  G  omemi  omema  metal  model  gain_model  solvent E0 rap host < ../data/input/nanosphere_eV.dat

if [ -n "$range_override" ]; then
    IFS=: read -r omi oma <<< "$range_override"
    if [ -z "$omi" ] || [ -z "$oma" ]; then
        echo "Error: -r|--range requires omemi:omema" >&2
        exit 1
    fi
    ex1=$(printf '%.6g' "$(echo "$omi - 0.3" | bc -l)")
    ex2=$(printf '%.6g' "$(echo "$oma + 0.3" | bc -l)")
elif [ "$metal" == "gold" ]; then
    ex1=1.6; ex2=2.8; omi=1.6; oma=2.8;
elif [ "$metal" == "silver" ]; then
    ex1=2.0; ex2=4.5; omi=2.2; oma=3.4;
fi

if [ -n "$rho_range_override" ]; then
    read -r -a rho_values <<< "$rho_range_override"
else
    rho_values=(0.4 0.5 0.6 0.7 0.8)
fi

if [ ${#rho_values[@]} -eq 0 ]; then
    echo "Error: -rho|--rho-range requires at least one value" >&2
    exit 1
fi

prev_rho=""
for rho in "${rho_values[@]}"; do
    if ! [[ "$rho" =~ ^[0-9]*\.?[0-9]+$ ]]; then
        echo "Error: -rho|--rho-range values must be positive numbers" >&2
        exit 1
    fi
    # ensure positive and increasing
    if (( $(echo "$rho <= 0" | bc -l) )); then
        echo "Error: -rho|--rho-range values must be positive" >&2
        exit 1
    fi
    if (( $(echo "$rho >= 1" | bc -l) )); then
        echo "Error: -rho|--rho-range values must be less than 1" >&2
        exit 1
    fi
    if [ -n "$prev_rho" ] && (( $(echo "$rho <= $prev_rho" | bc -l) )); then
        echo "Error: -rho|--rho-range values must be strictly increasing" >&2
        exit 1
    fi
    prev_rho="$rho"
done

echo "> Removing existing directory ../data/output/rho"
rm -fr ../data/output/rho
echo "> Creating directory ../data/output/rho"
mkdir -p ../data/output/rho

for rho in "${rho_values[@]}"; do
    echo "> rho = $rho"
    fro=(`./../bin/rho2ome_sp "$rho"`)
    omeG=${fro[0]}
    Gth=${fro[1]}
    GG=$(echo "1.2*$Gth" | bc -l)
    echo "> Setting the gain at the frohlich frequency ome_sp = $omeG eV"
    echo "> and setting the gain at 1.2*Gth = $GG ..."
    {
        echo "$a $Dome $omeG $GG $omi $oma $metal $model $gain_model $solvent $E0 $rho $host"
        echo "# a Dome omeG GG omi oma metal model gain_model solvent E0 rho host"
        echo "** WARNING: This file was automatically edited by $0."
        echo "** Do not edit this file manually until the process ends."
    } > ../data/input/nanosphere_eV.dat
    ./../bin/oap
    cp "../data/output/oGp/ome_p3.dat" "../data/output/rho/$rho.dat"
    echo "$omeG" > "../data/output/rho/omeB-$rho.dat"
    bash nano-shell-sketch.bash "$rho"
    mv "../img/nanoshell.png" "../data/output/rho/$rho.png"
done

echo "> Resetting the original input file..."
{
    echo "$a $Dome $ome21 $G $omemi $omema $metal $model $gain_model $solvent $E0 $rap $host"
    echo "# a Dome ome21 G omemi omema metal model gain_model solvent E0 rap host"
} > ../data/input/nanosphere_eV.dat
echo "> ...Done!"

rho_list="${rho_values[*]}"
xmax_list=()
ymax_list=()
for rho in "${rho_values[@]}"; do
    if [ ! -f "../data/output/rho/$rho.dat" ]; then
        echo "Warning: missing data file for rho=$rho, skipping stats for labels." >&2
        xmax_list+=("NaN")
        ymax_list+=("NaN")
        continue
    fi
    read x_max y_max < <(awk 'BEGIN{m=-1e30; x=0} {if($2>m){m=$2; x=$1}} END{printf "%.10g %.10g", x, m}' "../data/output/rho/$rho.dat")
    xmax_list+=("$x_max")
    ymax_list+=("$y_max")
done
rho_list_rev=()
xmax_rev=()
ymax_rev=()
for ((idx=${#rho_values[@]}-1; idx>=0; idx--)); do
    rho_list_rev+=("${rho_values[$idx]}")
    xmax_rev+=("${xmax_list[$idx]}")
    ymax_rev+=("${ymax_list[$idx]}")
done

xmax_list_str="${xmax_rev[*]}"
ymax_list_str="${ymax_rev[*]}"

dx_base=$(printf '%.10g' "$(echo "0.04*($oma-$omi)" | bc -l)")
declare -a dx_adjust
sketch_size=$dx_base
prev_edge=-1e30
prev_peak=""
# iterate from leftmost (lowest frequency) to rightmost
for idx in "${!xmax_rev[@]}"; do
    xval=${xmax_rev[$idx]}
    # start with default
    dx_val=$dx_base
    desired=$(echo "$xval + $dx_val" | bc -l)
    if [ -n "$prev_peak" ]; then
        gap=$(echo "$xval - $prev_peak" | bc -l)
        # if peaks are within 3*sketch size, give extra room
        if (( $(echo "$gap < $sketch_size*3" | bc -l) )); then
            dx_val=$(echo "$dx_base*2" | bc -l)
            desired=$(echo "$xval + $dx_val" | bc -l)
        fi
    fi
    # if overlapping previous sketch edge, double the offset once
    if (( $(echo "$desired <= $prev_edge" | bc -l) )); then
        dx_val=$(echo "$dx_base*2" | bc -l)
        desired=$(echo "$xval + $dx_val" | bc -l)
    fi
    dx_adjust[$idx]=$dx_val
    prev_edge=$(echo "$desired + $sketch_size" | bc -l)
    prev_peak=$xval
done

dx_list_str="${dx_adjust[*]}"

cat > "[nanoph-2024-0491]_rho_emi.gp" <<EOF
reset
omin=${omi}
omax=${oma}
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print \$col}' %s", row, col, file) )
r_list="${rho_list_rev[*]}"
n_r=words(r_list)
xmax_list="${xmax_list_str}"
ymax_list="${ymax_list_str}"
dx_list="${dx_list_str}"

omegaB(i) = at(sprintf("../data/output/rho/omeB-%s.dat", word(r_list,i)),1,1)
xmax(i) = real(word(xmax_list,i))
ymax_raw(i) = real(word(ymax_list,i))
dx(i) = real(word(dx_list,i))

# SETTING THE VISIBLE SPECTRUM IMAGE AT THE BOTTOM
set samples 200
set isosamples 2
k=omax-omin
set cbrange [omin:omax]
r(x)=x<2.1?1:x<2.4?-(x-2.4)/(2.4-2.1):x<2.8?0:x<3.3?(x-2.8)/(3.3-2.8):1
g(x)=x<1.9?0:x<2.1?(x-1.9)/(2.1-1.9):x<2.5?1:x<2.8?-(x-2.8)/(2.8-2.5):0
b(x)=x>2.5?1:x>2.4?-(x-2.4)/(2.4-2.5):0
f(x)=x<1.5?0:x<1.8?0.3+0.7*(x-1.6)/(1.8-1.6):x<3.?1:x<3.41?0.3+0.7*(3.3-x)/(3.3-3.):0
set palette functions f(k*gray+omin)*r(k*gray+omin),g(k*gray+omin),f(k*gray+omin)*b(k*gray+omin)
set palette functions f(k*gray+omin)*r(k*gray+omin),g(k*gray+omin),f(k*gray+omin)*b(k*gray+omin)
unset colorbox
# DONE

set term pdf color enhanced size 10cm, 8cm;
set output "../img/output/[nanoph-2024-0491]_rho.pdf"

set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03
set size 0.991,0.25
set pm3d map
unset ytics
set lmargin at screen 0.13
set rmargin at screen 0.97
set bmargin at screen 0.15
set tmargin at screen 0.2
set xlabel "ℏ{/Symbol w}_{em} (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

Isat=0.0926567

unset xlabel

set bmargin at screen 0.2
set tmargin at screen 0.94
unset xtics

set yrange [:0.2]
set ytics 0, 0.1

power=-1
div=10**power

scaleY(y)=y/(div*Isat)

# Determine dynamic y-range and place thumbnails/labels near peaks
max_plot_y = 0
do for [i=1:n_r] {
    py = scaleY(ymax_raw(i))
    if (py > max_plot_y) { max_plot_y = py }
}
set yrange [0:max_plot_y*1.1]

set for [i=1:n_r] pixmap (3+i) sprintf("../data/output/rho/%s.png", word(r_list,i)) at first (xmax(i)+dx(i)), first (0.9*scaleY(ymax_raw(i))) width screen 0.08
set for [i=1:n_r] label sprintf("{/Symbol r} = %s", word(r_list,i)) at first (xmax(i)+dx(i)), first (0.9*scaleY(ymax_raw(i))-0.01)

set ylabel "I_{em}/I_{sat}"
set xrange [omin:omax]
set label sprintf("x10^{%d}", power) at graph 0, 1.03

plot for [i=n_r:1:-1] sprintf("../data/output/rho/%s.dat", word(r_list,i)) u (\$1):(\$2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omegaB(i) t ""
unset multiplot
unset output
EOF

echo "> Producing the image..."
gnuplot "[nanoph-2024-0491]_rho_emi.gp"
echo "image ready in img/output/"
exit
