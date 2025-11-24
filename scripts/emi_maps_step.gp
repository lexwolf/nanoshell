#!/bin/bash
wG=abs(ARG1);
reset
set key bottom
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
file="../data/output/oGp/data4plot.dat" ;
# DONE
# row=1 ; col=1
# omin = at(file,row,col)
# row=1 ; col=2
# omax = at(file,row,col)
row=1 ; col=3
Gmax = at(file,row,col)
Gmax=abs(Gmax)

file="../data/input/nanosphere_eV.dat" ;
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)

set xlabel "Energy (eV)"
set ylabel "Quantity of gain G (n. u.)"
set xrange [omin:omax]
set yrange [0:Gmax]

set term pdf color enhanced size 11cm, 8cm;
set output "../img/output/oGp/map_isa_step.pdf"
set multiplot
set lmargin screen  0.15
set rmargin screen  0.75
set tmargin screen  0.95
set tmargin screen  0.95
set bmargin screen  0.15
set cblabel "|p_3|^2" offset 1.8,0.
set pm3d
unset surf
set view map
splot "../data/output/oGp/ome_G_p3.dat" u ($1):($2):($3) w l t ""
unset border 
unset tics
unset xlabel
unset ylabel
unset cblabel
set key textcolor variable
set arrow nohead from graph 0, first wG to graph 1, first wG lc rgb "red" dt 4 lw 2 front; 
plot "../data/output/oGp/iso_al.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "white" t "{/Symbol a}' < 0 and {/Symbol a}'' < 0", \
     0 w l lc rgb "black" t " "
unset multiplot
unset output
!convert -density 300 "../img/output/oGp/map_isa_step.pdf" "../img/output/oGp/map_isa_step.png"
!rm "../img/output/oGp/map_isa_step.pdf"

set border 
set tics
set xlabel
set ylabel
set cblabel

set output "../img/output/oGp/map_isk_step.pdf"
set multiplot
set lmargin screen  0.15
set rmargin screen  0.75
set tmargin screen  0.95
set tmargin screen  0.95
set bmargin screen  0.15
set cblabel "|p_3|^2" offset 1.8,0.
set pm3d
unset surf
set view map
splot "../data/output/oGp/ome_G_p3.dat" u ($1):($2):($3) w l t ""
# set border lc rgb "red"
unset border 
unset tics
unset xlabel
unset ylabel
unset cblabel
set key textcolor variable
set arrow nohead from graph 0, first wG to graph 0, first wG lc rgb "red" dt 2 lw 2;
plot "../data/output/oGp/iso_ka.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "gold" t "Re({/Symbol k}) > 0"
unset multiplot
unset output
!convert -density 300 "../img/output/oGp/map_isk_step.pdf" "../img/output/oGp/map_isk_step.png"
!rm "../img/output/oGp/map_isk_step.pdf"

