# Supplementary diagnostics: publication width, approximately 9-point labels.
# Run sm_figures.py to regenerate the plotting tables and full-sample metrics.
reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
load 'Pub_Ready/sm_transition_parameters.gp'
array gains[3] = [1.15, 1.20, 1.25]
array colors[4] = ['#1565C0', '#D97706', '#238567', '#8436A8']
array dashes[4] = [1, 2, 3, 4]
set encoding utf8
# Cairo applies a 2/3 text-size scale; 13.5-point input renders at 9 points.
set terminal pdfcairo enhanced color size 5.25in,6.75in font 'Tinos,13.5' linewidth 1.0
set border 15 linewidth 0.6
set xtics in mirror scale 0.55,0.3
set ytics in mirror scale 0.55,0.3
set mxtics 2
set mytics 2
set xrange [12:38]
set xtics 15,10,35
set samples 300
unset grid

do for [kind=1:2] {
    if (kind == 1) {
        set output '../LaTeX/sm_gain_switching_transition_mechanism.pdf'
    } else {
        set output '../LaTeX/sm_gain_switching_transition_cycles.pdf'
    }
    set multiplot
    # One compact shared legend above the six panels.
    set lmargin at screen 0.10
    set rmargin at screen 0.98
    set tmargin at screen 0.995
    set bmargin at screen 0.98
    unset border
    unset xtics
    unset ytics
    set yrange [0:1]
    set key at screen 0.54,0.991 center top horizontal samplen 1.7 spacing 1.0
    if (kind == 1) {
        set key maxcols 2 width 3
        plot 2 with lines lw 1.5 lc rgb '#8436A8' title '{/Tinos:Italic N}', \
             2 with lines lw 1.5 dt 2 lc rgb '#D97706' title '{/Tinos:Italic Ñ}', \
             2 with lines lw 0.9 dt 3 lc rgb '#666666' title '{/Tinos:Italic N}_{crit}', \
             2 with lines lw 0.8 dt 2 lc rgb '#AAAAAA' title 'Pump maximum'
    } else {
        set key maxcols 4 width 0
        plot for [cycle=1:4] 2 with lines lw 1.2 dt dashes[cycle] lc rgb colors[cycle] title sprintf('Cycle %d', cycle)
    }
    unset key
    set border 15 linewidth 0.6
    set xtics in mirror scale 0.55,0.3
    set ytics in mirror scale 0.55,0.3
    set xtics 15,10,35
    do for [row=1:3] {
        data = sprintf('Pub_Ready/sm_transition_%.2f_cycle4.dat', gains[row])
        set tmargin at screen (0.91-(row-1)*0.29)
        set bmargin at screen (0.705-(row-1)*0.29)
        set lmargin at screen 0.135
        set rmargin at screen 0.46
        set yrange [0.55:1.04]
        set ytics 0.6,0.2,1.0
        set format y '%.1f'
        set format x '%g'
        unset xlabel
        if (row == 3) { set xlabel 'Time within cycle (ps)' offset 0,0.4 }
        set ylabel 'Population inversion' offset 1.0,0
        set label 1 sprintf('(%s)',word('a c e',row)) at graph 0.04,0.91 front font 'Tinos:Bold,13.5'
        set label 2 sprintf('{/Tinos:Italic G}_{max}/{/Tinos:Italic G}_{th} = %.2f',gains[row]) at graph 0.5,1.13 center front
        set arrow 1 from graph 0,first 1/gains[row] to graph 1,first 1/gains[row] nohead back lw 0.9 dt 3 lc rgb '#666666'
        set arrow 2 from first 25,graph 0 to first 25,graph 1 nohead back lw 0.8 dt 2 lc rgb '#AAAAAA'
        if (kind == 1) {
            plot data using 1:3 with lines lw 1.5 dt 2 lc rgb '#D97706' notitle, \
                 data using 1:2 with lines lw 1.5 lc rgb '#8436A8' notitle
        } else {
            plot for [cycle=1:4] sprintf('Pub_Ready/sm_transition_%.2f_cycle%d.dat',gains[row],cycle) using 1:2 with lines lw 1.2 dt dashes[cycle] lc rgb colors[cycle] notitle
        }
        unset arrow 1
        unset label 2
        set lmargin at screen 0.635
        set rmargin at screen 0.96
        set yrange [0:12.5]
        set ytics 0,5,10
        set format y '%g'
        set ylabel '{/Tinos:Italic I}_{em}/{/Tinos:Italic I}_{sat} (×10^{-3})' offset 1.2,0
        set label 1 sprintf('(%s)',word('b d f',row)) at graph 0.04,0.91 front font 'Tinos:Bold,13.5'
        if (kind == 1) {
            set arrow 3 from first half_rises[row],first peaks[row]/2 to first half_falls[row],first peaks[row]/2 heads size screen 0.007,90 lw 1.2 lc rgb '#D97706' front
            plot data using 1:4 with filledcurves y=0 fs transparent solid 0.12 noborder lc rgb '#1565C0' notitle, \
                 data using 1:4 with lines lw 1.5 lc rgb '#1565C0' notitle, \
                 '+' using (peak_times[row]):(peaks[row]) every ::0::0 with points pt 7 ps 0.45 lc rgb '#1565C0' notitle
            unset arrow 3
        } else {
            plot for [cycle=1:4] sprintf('Pub_Ready/sm_transition_%.2f_cycle%d.dat',gains[row],cycle) using 1:4 with lines lw 1.2 dt dashes[cycle] lc rgb colors[cycle] notitle
        }
        unset arrow 2
    }
    unset label
    unset xlabel
    unset ylabel
    unset multiplot
    set output
}
