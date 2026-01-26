#!/bin/bash
set -o pipefail
export LC_NUMERIC="en_US.UTF-8"

# ---- float helpers (comparisons only) ----
f_le() { awk -v a="$1" -v b="$2" 'BEGIN{print (a<=b)?1:0}'; }
f_lt() { awk -v a="$1" -v b="$2" 'BEGIN{print (a<b)?1:0}'; }

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
    echo "                    (NOTE: quote the list, e.g. -rho \"0.4 0.6 0.8\")"
}

# ---- CLI parsing ----
while [ $# -gt 0 ]; do
    case "$1" in
        -c|--compile) compile_requested=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        -r|--range)
            if [ -z "${2:-}" ]; then
                echo "Error: -r|--range requires omemi:omema" >&2
                exit 1
            fi
            range_override="$2"
            shift 2
            ;;
        -rho|--rho-range)
            if [ -z "${2:-}" ]; then
                echo "Error: -rho|--rho-range requires a space-separated list (quote it)" >&2
                exit 1
            fi
            rho_range_override="$2"
            shift 2
            ;;
        --) shift; break ;;
        -*) echo "Unknown option: $1" >&2; show_help; exit 1 ;;
        *)  echo "Unexpected argument: $1" >&2; show_help; exit 1 ;;
    esac
done

if [ -n "$range_override" ] && [[ "$range_override" != *:* ]]; then
    echo "Error: -r|--range requires omemi:omema" >&2
    exit 1
fi

# ---- required files ----
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
    for item in "${missing[@]}"; do echo "  - $item"; done
    if [ "$compile_requested" != true ] && [ "$missing_exec" = true ]; then
        echo "> Hint: rerun with the -c flag to rebuild the executables."
    fi
    exit 1
fi

# ---- compile (optional) ----
if [ "$compile_requested" = true ]; then
    echo
    echo "> Compiling codes..."
    g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo || exit 1
    g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo || exit 1
    echo "> ...Done!"
    echo
fi

# ---- input file safety (ALWAYS restore) ----
input_file="../data/input/nanosphere_eV.dat"
[ -f "$input_file" ] || { echo "Error: missing $input_file" >&2; exit 1; }

tmp_backup="$(mktemp)"
cp -f "$input_file" "$tmp_backup" || exit 1
cleanup() {
    cp -f "$tmp_backup" "$input_file" >/dev/null 2>&1 || true
    rm -f "$tmp_backup" >/dev/null 2>&1 || true
}
trap cleanup EXIT

echo "> WARNING: $input_file will be temporarily overwritten during the calculations."

read a Dome ome21 G omemi omema metal model gain_model solvent E0 rap host < "$input_file"

# ---- range selection ----
if [ -n "$range_override" ]; then
    IFS=: read -r omi oma <<< "$range_override"
    if [ -z "$omi" ] || [ -z "$oma" ]; then
        echo "Error: -r|--range requires omemi:omema" >&2
        exit 1
    fi
    # arithmetic only
    ex1=$(printf '%.6g' "$(echo "$omi - 0.3" | bc -l)")
    ex2=$(printf '%.6g' "$(echo "$oma + 0.3" | bc -l)")
elif [ "$metal" == "gold" ]; then
    ex1=1.6; ex2=2.8; omi=1.6; oma=2.8
elif [ "$metal" == "silver" ]; then
    ex1=2.0; ex2=4.5; omi=2.2; oma=3.4
else
    echo "Error: metal='$metal' not handled and no -r override provided." >&2
    exit 1
fi

# ---- rho list ----
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
        echo "Error: rho values must be positive numbers" >&2
        exit 1
    fi
    (( $(f_le "$rho" 0) )) && { echo "Error: rho must be > 0" >&2; exit 1; }
    (( $(f_le 1 "$rho") )) && { echo "Error: rho must be < 1" >&2; exit 1; }
    if [ -n "$prev_rho" ] && (( $(f_le "$rho" "$prev_rho") )); then
        echo "Error: rho values must be strictly increasing" >&2
        exit 1
    fi
    prev_rho="$rho"
done

# ---- prepare output dir ----
echo "> Removing existing directory ../data/output/rho"
rm -fr ../data/output/rho
echo "> Creating directory ../data/output/rho"
mkdir -p ../data/output/rho

# ---- main loop ----
for rho in "${rho_values[@]}"; do
    echo "> rho = $rho"

    # safer than backticks/array splitting
    out="$(../bin/rho2ome_sp "$rho" 2>&1)"
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "Error: rho2ome_sp exited with code $rc for rho=$rho" >&2
        echo "Output was:" >&2
        echo "$out" >&2
        exit 1
    fi

    # take the first two numeric tokens from output (ignore colors / messages)
    read -r omeG Gth < <(awk '
        {
        for(i=1;i<=NF;i++){
            if($i ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/){
            n++; a[n]=$i;
            if(n==2){ print a[1],a[2]; exit }
            }
        }
        }
        END{ if(n<2) exit 1 }
    ' <<< "$out")

    if [ -z "${omeG:-}" ] || [ -z "${Gth:-}" ]; then
        echo "Error: rho2ome_sp did not provide two numeric values for rho=$rho" >&2
        echo "Raw output was:" >&2
        echo "$out" >&2
        exit 1
    fi


    GG=$(echo "1.2*$Gth" | bc -l)

    echo "> Setting the gain at the frohlich frequency ome_sp = $omeG eV"
    echo "> and setting the gain at 1.2*Gth = $GG ..."
    {
        echo "$a $Dome $omeG $GG $omi $oma $metal $model $gain_model $solvent $E0 $rho $host"
        echo "# a Dome omeG GG omi oma metal model gain_model solvent E0 rho host"
        echo "** WARNING: This file was automatically edited by $0."
        echo "** Do not edit this file manually until the process ends."
    } > "$input_file"

    ../bin/oap || { echo "Error: oap failed for rho=$rho" >&2; exit 1; }

    [ -f "../data/output/oGp/ome_p3.dat" ] || { echo "Error: missing ../data/output/oGp/ome_p3.dat after oap" >&2; exit 1; }
    cp "../data/output/oGp/ome_p3.dat" "../data/output/rho/$rho.dat" || exit 1
    echo "$omeG" > "../data/output/rho/omeB-$rho.dat"

    bash nano-shell-sketch.bash "$rho" || { echo "Error: nano-shell-sketch.bash failed for rho=$rho" >&2; exit 1; }
    [ -f "../img/nanoshell.png" ] || { echo "Error: nano-shell-sketch did not produce ../img/nanoshell.png" >&2; exit 1; }
    mv "../img/nanoshell.png" "../data/output/rho/$rho.png" || exit 1
done

echo "> Input file will now be restored (trap)."

# ---- find maxima for labels (numeric-only, initialized) ----
xmax_list=()
ymax_list=()

for rho in "${rho_values[@]}"; do
    dat="../data/output/rho/$rho.dat"
    if [ ! -f "$dat" ]; then
        echo "Warning: missing data file for rho=$rho, skipping stats for labels." >&2
        xmax_list+=("NaN")
        ymax_list+=("NaN")
        continue
    fi

    read x_max y_max < <(
      awk '
        BEGIN { m=-1e300; x="NaN"; seen=0 }
        $1 ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ &&
        $2 ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ {
            seen=1
            if($2>m){ m=$2; x=$1 }
        }
        END{
            if(!seen){ print "NaN NaN"; exit }
            printf "%.10g %.10g", x, m
        }' "$dat"
    )

    xmax_list+=("$x_max")
    ymax_list+=("$y_max")
done

# ---- reverse lists ----
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

# ---- compute dx adjustments (avoid NaN into bc) ----
dx_base=$(printf '%.10g' "$(echo "0.04*($oma-$omi)" | bc -l)")
declare -a dx_adjust
sketch_size=$dx_base
prev_edge=-1e30
prev_peak=""

for idx in "${!xmax_rev[@]}"; do
    xval=${xmax_rev[$idx]}

    if [[ "$xval" == "NaN" ]]; then
        dx_adjust[$idx]=$dx_base
        prev_peak=""
        continue
    fi

    dx_val=$dx_base
    desired=$(echo "$xval + $dx_val" | bc -l)

    if [ -n "$prev_peak" ]; then
        gap=$(echo "$xval - $prev_peak" | bc -l)
        if (( $(awk -v g="$gap" -v s="$sketch_size" 'BEGIN{print (g < 3*s)?1:0}') )); then
            dx_val=$(echo "$dx_base*2" | bc -l)
            desired=$(echo "$xval + $dx_val" | bc -l)
        fi
    fi

    if (( $(f_le "$desired" "$prev_edge") )); then
        dx_val=$(echo "$dx_base*2" | bc -l)
        desired=$(echo "$xval + $dx_val" | bc -l)
    fi

    dx_adjust[$idx]=$dx_val
    prev_edge=$(echo "$desired + $sketch_size" | bc -l)
    prev_peak=$xval
done

dx_list_str="${dx_adjust[*]}"

# ---- write gnuplot script (your original content preserved) ----
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
gnuplot "[nanoph-2024-0491]_rho_emi.gp" || { echo "Error: gnuplot failed" >&2; exit 1; }
echo "image ready in img/output/"
exit 0
