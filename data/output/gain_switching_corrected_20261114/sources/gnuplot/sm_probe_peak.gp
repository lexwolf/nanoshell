# Supplement: maximum dipolar response for the existing fixed-pump probe scan.
# Run: gnuplot Gain_Switching/gnuplot/sm_probe_peak.gp
reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] { if (script_path[i:i] eq '/') { last_slash = i } }
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
set encoding utf8
set terminal pdfcairo enhanced color size 5.25in,2.85in font 'Tinos,13.5' linewidth 1.1
set output '../LaTeX/sm_probe_peak_response.pdf'
set datafile commentschars '#'
set border 15 linewidth 1
set xtics in mirror
set ytics in mirror
set logscale xy 10
set xrange [1e-15:1]
set yrange [1e-3:1e4]
set xtics ('10^{-15}' 1e-15,'10^{-12}' 1e-12,'10^{-9}' 1e-9,'10^{-6}' 1e-6,'10^{-3}' 1e-3,'1' 1)
set ytics ('10^{-3}' 1e-3,'10^{-1}' 1e-1,'10^{1}' 10,'10^{3}' 1000)
set lmargin at screen 0.15
set rmargin at screen 0.97
set bmargin at screen 0.22
set tmargin at screen 0.96
set xlabel '{/Tinos:Italic E}_0 / {/Tinos:Italic E}_{sat}' offset 0,0.2
set ylabel '{/Tinos:Italic I}_{max} / {/Tinos:Italic I}_{sat}' offset 0.5,0
set key left top opaque samplen 2 spacing 1.1
array sources[3] = ['1.05Gth', '1.1Gth', '1.15Gth']
array gains[3] = ['1.05', '1.10', '1.15']
# Values read from each source header, not the instantaneous pump.
array esat[3] = [0.325413, 0.317931, 0.310942]
array colors[3] = ['#0072B2','#D55E00','#009E73']
plot for [i=1:3] sprintf('Data/ProbeSweep/middle/%s.dat',sources[i]) using 3:(($4/esat[i])**2) with lines lw 1.6 lc rgb colors[i] title sprintf('{/Tinos:Italic G}/{/Tinos:Italic G}_{th} = %s',gains[i])
set output
