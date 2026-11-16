# Rendering companion to scripts/gain_dependence.py. Run the Python entry
# point to recompute all metrics after replacing simulation files.
reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
data = 'Pub_Ready/gain_dependence_plot.dat'
set datafile missing 'nan'
stats data using 5 nooutput
has_flag = STATS_max > 0
array labels[3] = ['FWHM (ps)', '{/Tinos:Italic I}_{peak} / {/Tinos:Italic I}_{sat}', 'Peak time within cycle (ps)']
array colors[3] = ['#1565C0', '#D97706', '#8436A8']
set encoding utf8
do for [kind=1:2] {
    if (kind == 1) {
        set terminal pdfcairo enhanced color size 10.5in,4in font 'Tinos,11' linewidth 1.2
        set output 'Pub_Ready/gain_dependence.pdf'
    } else {
        set terminal pngcairo enhanced color size 3150,1200 font 'Tinos,33' linewidth 3.6
        set output 'Pub_Ready/gain_dependence.png'
    }
    point_factor = (kind == 1 ? 1.0 : 3.0)
    set multiplot
    set border 15 linewidth 0.8
    unset grid
    unset key
    set xtics in mirror scale 0.65,0.35
    set ytics in mirror scale 0.65,0.35
    set xrange [1.045:1.155]
    set xtics 1.05,0.02,1.15
    set format x '%.2f'
    set mxtics 2
    set xlabel '{/Tinos:Italic G}_{max} / {/Tinos:Italic G}_{th}'
    set label 90 'ω = ω_{th},  T = 50 ps,  E_0 / E_{sat} = 10^{-9}' at screen 0.52,0.965 center
    if (has_flag) { set label 92 'Draft: open markers denote unverified data' at screen 0.52,0.905 center textcolor rgb '#B34D00' }
    do for [panel=1:3] {
        set lmargin at screen (0.075+(panel-1)*0.325)
        set rmargin at screen (0.325+(panel-1)*0.325)
        set bmargin at screen 0.22
        set tmargin at screen 0.81
        set ylabel labels[panel] offset 0.5,0
        set label 1 sprintf('(%s)',word('a b c',panel)) at graph 0.06,0.92 front
        unset logscale y
        set autoscale y
        set ytics autofreq
        set format y '%g'
        set mytics 2
        if (panel == 2) {
            set logscale y
            set format y '10^{%L}'
            set mytics 10
        }
        plot data using 1:($5==0 ? column(panel+1) : 1/0) with linespoints pointtype 7 pointsize (0.65*point_factor) linewidth 1.5 linecolor rgb colors[panel] notitle, \
             data using 1:($5!=0 ? column(panel+1) : 1/0) with points pointtype 6 pointsize (1.1*point_factor) linewidth 1.6 linecolor rgb '#B34D00' notitle
        unset label 90
        unset label 92
    }
    unset multiplot
    set output
}
