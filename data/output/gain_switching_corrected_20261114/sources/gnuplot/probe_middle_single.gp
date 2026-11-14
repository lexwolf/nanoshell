# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Reproduce Images/probe_middle_single.pdf directly from the middle sweep.
# From Gain_Switching: gnuplot gnuplot/probe_middle_single.gp
# Paths are anchored to this script; invocation from other directories works.
# Source: ProbeSweep.ipynb cell 2 and DATA_IMAGE_MAP.md.
# Column 3 is E0/Esat (already normalized); column 1 is the time in ps
# required to reach 90% of the maximum dipole modulus.

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
set terminal pdfcairo enhanced color size 6.8052083in,4.2701667in font 'Tinos,18' linewidth 1.5
set output 'Pub_Ready/probe_middle_single.pdf'
set datafile commentschars '#'
unset grid
set border 15 linewidth 1.2
set xtics in mirror scale 0.65,0.35
set ytics in mirror scale 0.65,0.35
unset x2tics
unset y2tics
set logscale x 10
# Match the reference's 5% padding in log-x and linear-y coordinates.
set xrange [10**(-15.75):10**0.75]
set yrange [-1.54109365:34.89708065]
set xtics 1e-12,10000,1 format '10^{%L}'
set mxtics 10
set ytics 5
set mytics 5
set lmargin at screen 0.14
set rmargin at screen 0.985
set bmargin at screen 0.19
set tmargin at screen 0.975
set xlabel '{/Tinos:Italic E}_0 / {/Tinos:Italic E}_{sat}' font 'Tinos,18' offset 0,0.1
set ylabel '{/Tinos:Italic t}_{on} (ps)' font 'Tinos,18' offset 0.5,0
set style line 10 linecolor rgb '#BBBBBB' linewidth 0.8
set key right top reverse Left box linestyle 10 opaque samplen 2.5 spacing 1.2

array gains[3] = ['1.05Gth', '1.1Gth', '1.15Gth']
array titles[3] = ['1.05', '1.10', '1.15']
array colors[3] = ['#0072B2', '#D55E00', '#009E73']
plot for [p=1:3] sprintf('Data/ProbeSweep/middle/%s.dat', gains[p]) \
    using 3:1 with lines linewidth 1.8 linecolor rgb colors[p] \
    title sprintf('%s {/Tinos:Italic G}_{th}', titles[p])
set output
