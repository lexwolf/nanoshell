reset
set terminal pdfcairo enhanced color size 28cm,9cm font "Arial,24"

repo = (system("test -d scripts -a -d data -a -d img; echo $?") eq "0") ? "." : ".."
data = repo."/data/output/intensity_vs_gain/intensity_vs_gain.dat"

set border lw 1.6
set tics out nomirror
set grid lc rgb "#d0d0d0" lw 0.6
set key left top spacing 1.05 samplen 2.2

set xlabel "Gain G/G_{th}"
set ylabel "Emission intensity / I_{sat}"
set logscale y
set xrange [0:2]
set yrange [:100]
set xtics 0, 0.2
set ytics 1e-14, 1.e4
set format y "10^{%T}"

set style line 1 lc rgb "#e9551a" lw 6
set style line 2 lc rgb "#2a9d55" lw 6
set style line 3 lc rgb "#7b3294" lw 8 dt 2
set style line 4 lc rgb "#666666" lw 1.5 dt 2

set obj rect from graph 0, graph 0 to graph 1, graph 1
set arrow 1 from 1, graph 0 to 1, graph 1 nohead ls 4 front

set output repo."/img/output/intensity_vs_gain0.pdf"
plot data using 2:3 with lines ls 1 title "|{/Symbol a}_{SS}E_0|^2/I_{sat}"
unset output

set output repo."/img/output/intensity_vs_gain1.pdf"
plot data using 2:3 with lines ls 1 title "|{/Symbol a}_{SS}E_0|^2/I_{sat}", \
     data using 2:4 with lines ls 2 title "|p_{num}|^2/I_{sat}"
unset output

set output repo."/img/output/intensity_vs_gain.pdf"
plot data using 2:3 with lines ls 1 title "|{/Symbol a}_{SS}E_0|^2/I_{sat}", \
     data using 2:4 with lines ls 2 title "|p_{num}|^2/I_{sat}", \
     data using 2:5 with lines ls 3 title "|p_{analytic}|^2/I_{sat}"
unset output
