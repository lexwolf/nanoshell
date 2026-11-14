# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Reproduce Images/freq_sweep_LR_panels.png as a vector PDF.
# Run from Gain_Switching: gnuplot gnuplot/freq_sweep_LR_panels.gp
# Paths are anchored to this script, independent of the working directory.
# Source: DATA_IMAGE_MAP.md and FreqSweep.ipynb cell 4.
# Column 5 = emission energy (eV); column 1 = 90%-dipole turn-on time (ps).
# Left/Right identify study-frequency branches. Preserve each file's row order.

reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
system 'mkdir -p Pub_Ready'

set encoding utf8
set terminal pdfcairo enhanced color size 10.1166667in,5.1in font 'Tinos,26.4' linewidth 1.5
set output 'Pub_Ready/freq_sweep_LR_panels.pdf'
set datafile commentschars '#'
unset grid
set border 15 linewidth 1.2
set xtics in mirror scale 0.65,0.35
set ytics in mirror scale 0.65,0.35
unset x2tics
unset y2tics
set ytics 20
set mxtics 5
set mytics 5
set yrange [0:100]
set xlabel 'ℏ{/Tinos:Italic ω}_{em} (eV)' font 'Tinos,26.4' offset 0,0.2
set style line 10 linecolor rgb '#BBBBBB' linewidth 0.8
set key left top reverse Left box linestyle 10 opaque samplen 1 spacing 1.1

array gains[3] = ['1.05Gth', '1.1Gth', '1.15Gth']
array titles[3] = ['1.05', '1.10', '1.15']
array colors[3] = ['#0072B2', '#D55E00', '#009E73']
array letters[3] = ['(a)', '(b)', '(c)']
# Complete branch extents plus the reference's 5% horizontal padding.
array xmin[3] = [2.8122, 2.8122, 2.8122]
array xmax[3] = [2.82772, 2.8346, 2.83987]

set multiplot
do for [p=1:3] {
    left = 0.105 + (p-1)*0.300
    set lmargin at screen left
    set rmargin at screen left+0.265
    set bmargin at screen 0.23
    set tmargin at screen 0.85
    pad = 0.05*(xmax[p]-xmin[p])
    set xrange [xmin[p]-pad:xmax[p]+pad]
    if (p < 3) {
        set xtics 2.815,(p == 1 ? 0.005 : 0.01),2.835 format '%.3f'
    } else {
        set xtics 0.01 format '%.2f'
    }
    set title sprintf('{/Tinos:Italic G} = %s {/Tinos:Italic G}_{th}', titles[p]) font 'Tinos,26.4' offset 0,0.1
    if (p == 1) {
        set format y '%g'
        set ylabel '{/Tinos:Italic t}_{on} (ps)' font 'Tinos,26.4' offset 0.5,0
    } else {
        set format y ''
        unset ylabel
    }
    set label 1 letters[p] at graph 0.92,0.10 right front font 'Tinos:Bold,26.4'
    plot sprintf('Data/FreqSweep/%s/Left.dat', gains[p]) using 5:1 with lines linewidth 1.8 dashtype 1 linecolor rgb colors[p] title 'Left', \
         sprintf('Data/FreqSweep/%s/Right.dat', gains[p]) using 5:1 with lines linewidth 1.8 dashtype 2 linecolor rgb colors[p] title 'Right'
}
unset multiplot
set output
