# Active OMEx plot: original publication styling; inputs are NEW corrected calculations.
# Source reference retained in scripts/[OMEx]_reference/gnuplot/.
# Normalization comes from validated run metadata; intensity limits adjusted to prevent clipping.
# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Reproduce Images/freq_sweep_panels.pdf from the original numeric data.
# Run with gnuplot, e.g. from Gain_Switching:
#   gnuplot gnuplot/freq_sweep_panels.gp
# Paths are anchored to this script, so invocation from another directory works.
# Source: DATA_IMAGE_MAP.md and FreqSweep.ipynb cell 0.
# Column 2 = study energy (eV); column 1 = 90%-dipole turn-on time (ps).
# Keep original row order and clip the ordinate, as in the reference figure.

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
system 'mkdir -p img/output/[OMEx]_gain_switching'

set encoding utf8
set terminal pdfcairo enhanced color size 10.1166667in,5.1in font 'Tinos,26.4' linewidth 1.5
set output 'img/output/[OMEx]_gain_switching/[OMEx]_freq_sweep_panels.pdf'
set datafile commentschars '#'
unset key
unset grid
set border 15 linewidth 1.2
set xtics in mirror scale 0.65,0.35
set ytics in mirror scale 0.65,0.35
unset x2tics
unset y2tics
set xtics 0.2 format '%.1f'
set ytics 20
set mxtics 5
set mytics 5
set yrange [0:100]
set xlabel 'ℏ{/Tinos:Italic ω} (eV)' font 'Tinos,26.4' offset 0,0.2

array gains[3] = ['1.05Gth', '1.1Gth', '1.15Gth']
array titles[3] = ['1.05', '1.10', '1.15']
array colors[3] = ['#0072B2', '#D55E00', '#009E73']
array letters[3] = ['(a)', '(b)', '(c)']
# Matplotlib's reference uses a 5% x margin around the complete data extent.
array xmin[3] = [2.51003, 2.45272, 2.41524]
array xmax[3] = [3.1454, 3.21648, 3.26449]

set multiplot
do for [p=1:3] {
    left = 0.105 + (p-1)*0.300
    set lmargin at screen left
    set rmargin at screen left+0.265
    set bmargin at screen 0.23
    set tmargin at screen 0.85
    pad = 0.05*(xmax[p]-xmin[p])
    set xrange [xmin[p]-pad:xmax[p]+pad]
    set title sprintf('{/Tinos:Italic G} = %s {/Tinos:Italic G}_{th}', titles[p]) font 'Tinos,26.4' offset 0,0.1
    if (p == 1) {
        set format y '%g'
        set ylabel '{/Tinos:Italic t}_{on} (ps)' font 'Tinos,26.4' offset 0.5,0
    } else {
        set format y ''
        unset ylabel
    }
    set label 1 letters[p] at graph 0.92,0.10 right front font 'Tinos:Bold,26.4'
    plot sprintf('data/output/gain_switching_corrected_20261114/publication_plot_data/Data/FreqSweep/%s/%s.dat', gains[p], gains[p]) using 2:1 with lines linewidth 1.8 linecolor rgb colors[p]
}
unset multiplot
set output
