# Active OMEx plot: original publication styling; inputs are NEW corrected calculations.
# Source reference retained in scripts/[OMEx]_reference/gnuplot/.
# Normalization comes from validated run metadata; intensity limits adjusted to prevent clipping.
# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Linear-scale comparison variant. Uses the existing validated metrics.
# Run with gnuplot; writes separate outputs and leaves the manuscript unchanged.
reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
load 'scripts/[OMEx]_publication_parameters.gp'
load 'scripts/[OMEx]_gain_dependence_fits.gp'
data = 'data/output/gain_switching_corrected_20261114/publication_plot_data/Pub_Ready/gain_dependence_linear_extended_plot.dat'
set datafile missing 'nan'
stats data using 5 nooutput
has_flag = STATS_max > 0
array labels[3] = ['FWHM (ps)', '{/Tinos:Italic I}_{peak} / {/Tinos:Italic I}_{sat}', 'Peak time (ps)']
array colors[3] = ['#1565C0', '#D97706', '#8436A8']
set encoding utf8
do for [kind=1:2] {
    if (kind == 1) {
        set terminal pdfcairo enhanced color size 10.5in,6in font 'Tinos,27' linewidth 1.2
        set output 'img/output/[OMEx]_gain_switching/[OMEx]_gain_dependence_linear_extended.pdf'
    } else {
        set terminal pngcairo enhanced color size 3150,1800 font 'Tinos,81' linewidth 3.6
        set output 'img/output/[OMEx]_gain_switching/[OMEx]_gain_dependence_linear_extended.png'
    }
    point_factor = (kind == 1 ? 1.0 : 3.0)
    set multiplot
    set border 15 linewidth 0.8
    unset grid
    unset key
    set xtics in mirror scale 0.65,0.35
    set ytics in mirror scale 0.65,0.35
    set xrange [1.04:1.51]
    set xtics 1.10,0.20,1.50
    set format x '%.2f'
    set mxtics 2
    set xlabel '{/Tinos:Italic G}_{max} / {/Tinos:Italic G}_{th}'
    set label 90 'ω = ω_{th},  T = 50 ps,  E_0 / E_{sat} = 10^{-9}' at screen 0.52,0.95 center
    if (has_flag) { set label 92 'Draft: open markers denote unverified data' at screen 0.52,0.905 center textcolor rgb '#B34D00' }
    do for [panel=1:3] {
        set lmargin at screen (0.12+(panel-1)*0.32)
        set rmargin at screen (0.345+(panel-1)*0.32)
        set bmargin at screen 0.24
        set tmargin at screen 0.80
        set ylabel labels[panel] offset 0.5,0
        set label 1 sprintf('(%s)',word('a b c',panel)) at graph (panel == 3 ? 0.83 : 0.06),0.92 front
        unset logscale y
        set autoscale y
        set ytics autofreq
        set format y '%g'
        set mytics 2
        if (panel == 2) {
            set yrange [0:0.22]
            set format y '%.2f'
            set ytics 0.10
        }
        if (panel == 1) { set ytics 4 }
        if (panel == 3) { set ytics 2 }
        plot data using 1:($5==0 ? column(panel+1) : 1/0) with linespoints pointtype 7 pointsize (0.65*point_factor) linewidth 1.5 linecolor rgb colors[panel] notitle, \
             data using 1:($5!=0 ? column(panel+1) : 1/0) with points pointtype 6 pointsize (1.1*point_factor) linewidth 1.6 linecolor rgb '#B34D00' notitle
        if (panel == 2) {
#             set label 1 '(b)' at screen 0.4535,0.7552 front
            unset label 1
            unset xlabel
            unset ylabel
#             set object 99 rect from screen 0.46,0.31 to screen 0.625,0.47 behind fc rgb 'white' fs solid 1 noborder
            set lmargin at screen 0.555
            set rmargin at screen 0.66
            set bmargin at screen 0.26
            set tmargin at screen 0.445
            set logscale y
            set xrange [1.045:1.185]
            set yrange [1e-8:0.12]
            unset xtics
#             set xtics ('1.05' 1.05,'1.15' 1.15) font (kind == 1 ? 'Tinos,18' : 'Tinos,54') offset 0,0.7
            set x2tics mirror ('1.05' 1.05,'1.15' 1.15) font (kind == 1 ? 'Tinos,18' : 'Tinos,54') offset 0,-0.7
            set ytics ('10^{-8}' 1e-8,'10^{-6}' 1e-6,'10^{-4}' 1e-4,'10^{-2}' 1e-2) font (kind == 1 ? 'Tinos,18' : 'Tinos,54')
            set format y '%g'
            set mxtics 1
            set mytics 2
            set label 93 sprintf('m ≈ %.1f',omex_low_m) at graph 0.7,0.12 right front textcolor rgb '#1B5E20' font (kind == 1 ? 'Tinos,18' : 'Tinos,54')
            plot omex_low_fit(x) with lines dt 1 lw 3 lc rgb '#1B5E20' notitle, \
                 data using 1:3 with linespoints pt 7 ps (0.65*point_factor) lw 1 lc rgb colors[2] notitle
            unset label 93
            unset x2tics
            unset object 99
            unset label 1
            unset logscale y
            set xrange [1.04:1.51]
            set xtics 1.10,0.20,1.50 font '' offset 0,0
            set ytics font ''
            set mxtics 2
            set xlabel '{/Tinos:Italic G}_{max} / {/Tinos:Italic G}_{th}'
        }
        unset label 90
        unset label 92
    }
    unset multiplot
    set output
}
