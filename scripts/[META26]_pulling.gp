reset
omin=1.8
omax=3.65
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
omegaB = at("../data/output/pulling/omeB.dat",1,1)
rho = at("../data/input/nanosphere_eV.dat",1,12)

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

set term pdf color enhanced size 12cm, 10cm font "Arial,18";
set output "../img/output/[META26]_pulling.pdf"

set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03
set size 0.991,0.25
set pm3d map
unset ytics
set lmargin at screen 0.18
set rmargin at screen 0.98
set bmargin at screen 0.15
set tmargin at screen 0.2
set xlabel "ℏ{/Symbol w}_{em} (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

Isat=0.0926567

unset xlabel

set bmargin at screen 0.2
set tmargin at screen 0.94
unset xtics

set yrange [:1]
set ytics 0, 0.2

power=-2
div=10**power

scaleY(y)=y/(div*Isat)

set pixmap 3 sprintf("../data/output/pulling/%s.png", rho) at graph 0.21, graph 0.8 width screen 0.09
set label sprintf("{/Symbol r} = %s", rho) at graph 0.05, graph 0.875 left

set ylabel "I_{em}/I_{sat}"
set xrange [omin:omax]
set label sprintf("x10^{%d}", power) at graph 0, 1.03

plot "../data/output/pulling/ome_p3.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omegaB t "", \
     "../data/output/ns_intensity_SS.dat" u ($1):($5/(div*Isat)) w l lw 2 lc palette cb omegaB t "";
unset multiplot
unset output
