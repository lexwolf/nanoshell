reset
omin=1.8
omax=3.65
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
# DONE

file="../data/output/rho/omeB-0.4.dat" ;
row=1 ; col=1
omeB04 = at(file,row,col)

file="../data/output/rho/omeB-0.5.dat" ;
row=1 ; col=1
omeB05 = at(file,row,col)

file="../data/output/rho/omeB-0.6.dat" ;
row=1 ; col=1
omeB06 = at(file,row,col)

file="../data/output/rho/omeB-0.7.dat" ;
row=1 ; col=1
omeB07 = at(file,row,col)

file="../data/output/rho/omeB-0.8.dat" ;
row=1 ; col=1
omeB08 = at(file,row,col)

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
# 
set output "../img/output/[nanoph-2024-0491]_rho.pdf"

set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03
set size 0.991,0.25
# set border lt 1 lc 0 lw 3
set pm3d map
unset ytics
set lmargin at screen 0.13
set rmargin at screen 0.99
set bmargin at screen 0.15
set tmargin at screen 0.2
set xlabel "ℏ{/Symbol w}_{em} (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

Isat=0.0926567

xpos=3.16
ypos=0.04
xlab=xpos+0.22
ylab=ypos+0.09
set pixmap 4 "../data/output/rho/0.4.png" 
set pixmap 4 at first xpos, graph ypos width screen 0.09
set label at first xlab, graph ylab "{/Symbol r} = 0.4"

xpos=2.99
ypos=0.22
xlab=xpos+0.22
ylab=ypos+0.09
set pixmap 5 "../data/output/rho/0.5.png" 
set pixmap 5 at first xpos, graph ypos width screen 0.09
set label at first xlab, graph ylab "{/Symbol r} = 0.5"


xpos=2.8
ypos=0.42
xlab=xpos+0.22
ylab=ypos+0.09
set pixmap 6 "../data/output/rho/0.6.png" 
set pixmap 6 at first xpos, graph ypos width screen 0.09
set label at first xlab, graph ylab "{/Symbol r} = 0.6"

xpos=2.6
ypos=0.6
xlab=xpos+0.22
ylab=ypos+0.09
set pixmap 7 "../data/output/rho/0.7.png" 
set pixmap 7 at first xpos, graph ypos width screen 0.09
set label at first xlab, graph ylab "{/Symbol r} = 0.7"

xpos=2.27
ypos=0.82
xlab=xpos+0.22
ylab=ypos+0.09
set pixmap 8 "../data/output/rho/0.8.png" 
set pixmap 8 at first xpos, graph ypos width screen 0.09
set label at first xlab, graph ylab "{/Symbol r} = 0.8"

unset xlabel

set bmargin at screen 0.2
set tmargin at screen 0.94
unset xtics

set yrange [:0.2]
set ytics 0, 0.1

power=-1
div=10**power

set ylabel "I_{em}/I_{sat}"
set xrange [omin:omax]
set label sprintf("x10^{%d}", power) at graph 0, 1.03

plot "../data/output/rho/0.8.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB08 t "",\
     "../data/output/rho/0.7.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB07 t "",\
     "../data/output/rho/0.6.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB06 t "",\
     "../data/output/rho/0.5.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB05 t "",\
     "../data/output/rho/0.4.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB04 t "";
unset multiplot
unset output
