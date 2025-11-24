reset
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
file="../data/input/nanosphere_eV.dat" ;
# DONE

# READING SPECTRUM RANGE (omin omax) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)
# DONE

# READING GAIN QUANTITY (G) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=4
G = abs(at(file,row,col))
# DONE

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

set term pdf color enhanced size 10cm, 8cm;
set output "../img/output/spectrum.pdf"

# set term postscript eps enhanced color
sz=1.
set size sz
# set output "../img/output/spectrum.eps"
set multiplot

# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03*sz
set size 0.991*sz,0.25*sz
# set border lt 1 lc 0 lw 3
set pm3d map
set lmargin screen 0.2*sz
set rmargin screen 0.88*sz
unset ytics
set xlabel "Energy (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

# PLOT SETTINGS
set y2label "Polarizability"
set border lt 1 lc 0 lw 1
unset xlabel
set origin 0,0.171*sz
set size sz,0.85*sz
unset xtics
set xzeroaxis lt 1 lc 0

set yrange[-0.00005:0.00014]
set yrange[-0.00003:0.00006]

set ytics nomirror;
set y2tics;
set link y2 via y*1.e5 inverse y*1.e-5;
set ylabel "|q_0|^2 (n. u.)"
set y2label "Polarizability"


plot[omin:omax]  "../data/output/stationary.dat" u 1:3 w l lw 2 lc rgb "red" axis x1y2 t "" , \
                 "../data/output/stationary.dat" w l lw 2 lc rgb "black" axis x1y2 t "", \
                 "../data/output/ns_funct_ints_SS.dat" u 1:3 w l lw 2 lc 1 t "non-shifted |q_0|^2", \
                 "../data/output/right.dat" w l lc 2 lw 2 t "shifted |q_0|^2", \
                 NaN w l lw 2 lc rgb "red" axis x1y2 t "Im[{/Symbol a}]", NaN w l lw 2 lc rgb "black" axis x1y2 t "Re[{/Symbol a}]";
