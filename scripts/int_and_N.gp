reset
set term pdf color enhanced size 20cm, 10cm font ",22";
set output '../img/output/int_and_N.pdf';
set style line 6 lt 1 lc 0 lw 3;
set style line 7 lt 1 lc 8 lw 3;
set style line 8 lc rgb "purple" lw 4;
set style line 9 lc rgb "black" lw 4;
set style line 10 lt 1.5 lc rgb "magenta" lw 3;
set style line 11 lt 1.5 lc rgb "dark-green" lw 2;

set style fill transparent solid 0.4

power=-5
div=10**power

set xrange [0:600]
# set ytics -0.01, 0.01, 0.01
set key left
# set size 1,0.319

# set bmargin 4
# set lmargin 5.5
# set rmargin 7

set tmargin 1.5

set ytics nomirror
set ylabel "Field Intensity in the Core"
set y2tics tc rgb "magenta"
set y2label "Population Inversion" tc rgb "magenta"
set label sprintf("x10^{%d}", power) at graph 0, 1.05

# set y2range [98:100]
set y2tics 0.98, 0.01
set yrange [0:14]
set y2range [.98:1]
set format y2 "%.2f"
set xlabel "Time (ps)"
plot "../data/output/numtime.dat" u ($1):(($7*$7+$8*$8)/div) w filledcurve x1 lt 1 lw 2 lc rgb "dark-green" t "|p_0|^2",\
     "../data/output/numtime.dat" u ($1):(($7*$7+$8*$8)/div) w l lt 1 lw 2 lc rgb "dark-green" t "",\
     "../data/output/numtime.dat" u ($1):($6) axis x1y2 w l ls 10 t "N";
