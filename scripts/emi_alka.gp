reset
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
file="../data/input/nanosphere_eV.dat" ;
# DONE

# READING SPECTRUM RANGE (omin omax) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=3
omeG = at(file,row,col)
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)

p3mx=0.001
cnvrsn0=p3mx/20.
# DONE

# READING GAIN QUANTITY (G) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=4
G = abs(at(file,row,col))
# DONE

file="../data/output/oGp/ome_ex.dat" ;
row=1 ; col=1
omel = at(file,row,col)
row=1 ; col=2
omer =at(file,row,col)
row=1 ; col=3
omeB =at(file,row,col)
row=1 ; col=4
kexl = at(file,row,col)
row=1 ; col=5
kexr =at(file,row,col)

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

sz=1.
set size sz

set output "../img/output/emi_alka.pdf"

x1=0.12;  y1=0.95;
x0=0.06;  y0=0.95;
if (abs(omax - omeG) < abs(omin-omeG)) \
    set key left bottom; \
    else  x1=0.65;  y1=0.05;\
          x0=0.59;  y0=0.05;
    
set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03*sz
set size 0.991*sz,0.26*sz
# set border lt 1 lc 0 lw 3
set pm3d map
set lmargin screen 0.15*sz
set rmargin screen 0.8*sz
unset ytics
set xlabel "Energy (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

# PLOT SETTINGS
set ytics nomirror
set y2tics
set border lt 1 lc 0 lw 1
unset xlabel
set origin 0,0.171*sz
set size sz,0.85*sz
unset xtics
set xzeroaxis lt 1 lc 0
# set key over
set ylabel "Normalized Steady State Polarizability"
set y2label "Amplification Eigenvalue"
set yrange [-15:25]
set link y2 via cnvrsn0*y inverse y/cnvrsn0
set style rect fc lt -1 fs solid 0.3 noborder
# set obj rect from first omel, graph 0 to first omer, graph 1 fc rgb "black" 
set style rect fc lt -1 fs solid 0.3 border
# set obj rect at graph x0, graph y0 size graph 0.1, graph 0.05 fc rgb "black" 
# set label at graph x1, graph y1 "{/Symbol a}' < 0 and {/Symbol a}'' < 0"
plot[omin:omax]  "../data/output/oGp/ome_al.dat" u ($1):($4) axis x1y2  w l lw 2 lc rgb "orange" t "Re({/Symbol k})",\
                 ""  u ($1):($5/1000) axis x1y2  w l lw 2 lc rgb "purple" t "Im({/Symbol k})×10^{-3}",\
                 "" u 1:3 w l lw 2 lc rgb "red" t "{/Symbol a}''/4{/Symbol p}a^3";

unset multiplot
unset output
