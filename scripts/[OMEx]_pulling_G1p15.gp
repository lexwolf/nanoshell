# Active OMEx plot: original publication styling; inputs are NEW corrected calculations.
# Source reference retained in scripts/[OMEx]_reference/gnuplot/.
# Normalization comes from validated run metadata; intensity limits adjusted to prevent clipping.
# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
reset
# Resolve project-relative data/output paths independently of the launch directory.
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
load 'data/output/pulling_G1p15_20261120/plot_data/parameters.gp'
system 'mkdir -p img/output/[OMEx]_gain_switching'
omin=1.8
omax=3.6
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
omegaB = at("data/output/pulling_G1p15_20261120/plot_data/omeB.dat",1,1)
rho = at("data/output/pulling_G1p15_20261120/plot_data/nanosphere_eV.dat",1,12)

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

# Find the two envelope intersections at each illustrated intensity level.
# Levels are in the displayed ordinate (Iem/Isat multiplied by 10^2).
array arrow_levels[3] = [1.0804516033650513, 3.5972682794389352, 5.6183483374982659]
array arrow_left[3]
array arrow_right[3]
array arrow_target[3]
Isat=candidate_Isat
power=-2
div=10**power
previous_x=0.0
previous_y=0.0
have_previous=0
intersection(x,y) = (hit=(have_previous && (previous_y-level)*(y-level)<0 ? previous_x+(level-previous_y)*(x-previous_x)/(y-previous_y) : 1/0), previous_x=x, previous_y=y, have_previous=1, hit)
do for [j=1:3] {
    level=arrow_levels[j]
    have_previous=0
    stats 'data/output/pulling_G1p15_20261120/plot_data/ns_intensity_SS.dat' using (intersection($1,$5/(div*Isat))) nooutput
    arrow_left[j]=STATS_min
    arrow_right[j]=STATS_max
    # The rightmost intersection is on the emission curve; the left edge
    # of the filled band is an artificial closure near the threshold energy.
    have_previous=0
    stats 'data/output/pulling_G1p15_20261120/plot_data/ome_p3.dat' using (intersection($1,$2/(div*Isat))) nooutput
    arrow_target[j]=STATS_max
    print sprintf('Arrow level %.2f: %.6f -> %.6f <- %.6f eV', level, arrow_left[j], arrow_target[j], arrow_right[j])
}

set encoding utf8
set terminal pdfcairo enhanced color size 7in,4.5in font 'Tinos,18' linewidth 1.5
set border 15 linewidth 1
unset grid
set xtics in mirror scale 0.65,0.35
set ytics in mirror scale 0.65,0.35
set output 'img/output/[OMEx]_gain_switching/[OMEx]_pulling_G1p15.pdf'

set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03
set size 0.991,0.25
set pm3d map
unset ytics
set lmargin at screen 0.16
set rmargin at screen 0.975
set bmargin at screen 0.19
set tmargin at screen 0.25
set xlabel "Energy (eV)" font "Tinos,18" offset 0,0
set xtics 0.2 offset 0,-0.35
set mxtics 2
splot[omin:omax] x t ""
# DONE

Isat=candidate_Isat

unset xlabel

set bmargin at screen 0.25
set tmargin at screen 0.91
unset xtics

set yrange [0:7]
set ytics 0, 2
set mytics 2
set format y '%.1f'

power=-2
div=10**power

scaleY(y)=y/(div*Isat)


set ylabel "{/Tinos:Italic I}_{em} / {/Tinos:Italic I}_{sat}" font "Tinos,18" offset 0.5,0
set xrange [omin:omax]
set label sprintf("×10^{%d}", power) at graph 0, 1.03 font "Tinos,18"

# Palette-colored guides and gradient arrows, pointing to the emission curve.
# Each shaft segment samples the same palette coordinate as the spectrum strip.
segments=120
arrow_id=100
do for [j=1:3] {
    level=arrow_levels[j]
    target_x=arrow_target[j]
    do for [side=1:2] {
        start_x=(side==1 ? arrow_left[j] : arrow_right[j])
        arrow_id=arrow_id+1
        set arrow arrow_id from first start_x,0 to first start_x,level nohead front lw 1.6 lc palette cb start_x
        do for [seg=0:segments-1] {
            x0=start_x+(target_x-start_x)*seg/segments
            x1=start_x+(target_x-start_x)*(seg+1)/segments
            arrow_id=arrow_id+1
            set arrow arrow_id from first x0,level to first x1,level nohead front lw 1.6 lc palette cb ((x0+x1)/2.0)
        }
        arrow_id=arrow_id+1
        head_start=target_x+(side==1 ? -0.07 : 0.07)
        set arrow arrow_id from first head_start,level to first target_x,level head filled size screen 0.018,20 fixed front lw 1.6 lc palette cb target_x
    }
}

plot "data/output/pulling_G1p15_20261120/plot_data/ome_p3.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omegaB t "", \
     "data/output/pulling_G1p15_20261120/plot_data/ns_intensity_SS.dat" u ($1):($5/(div*Isat)) w l lw 2 lc palette cb omegaB t "";
unset multiplot
unset output
