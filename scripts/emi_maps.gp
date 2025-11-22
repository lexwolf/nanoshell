#!/bin/bash
reset
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
file="../data/output/oGp/data4plot.dat" ;
# DONE
row=1 ; col=1
omin = at(file,row,col)
row=1 ; col=2
omax = at(file,row,col)
row=1 ; col=3
Gmax = at(file,row,col)
Gmax=abs(Gmax)

set xlabel "Energy (eV)"
set ylabel "Quantity of gain G (n. u.)"
set xrange [omin:omax]
set yrange [0:Gmax]

set key bottom
set term pdf color enhanced size 11cm, 8cm;
set lmargin screen  0.15
set rmargin screen  0.75
set tmargin screen  0.95
set tmargin screen  0.95
set bmargin screen  0.15
set output "../img/output/oGp/iso_alka.pdf"
plot "../data/output/oGp/iso_al.dat" u ($1):($2) w filledcurves fs transparent solid 0.3 lw 2 lc rgb "black" t "{/Symbol a}' < 0 and {/Symbol a}'' < 0", \
     "../data/output/oGp/iso_ka.dat" u ($1):($2) w filledcurves fs transparent solid 0.3 lw 2 lc rgb "gold" t "Re({/Symbol k}) > 0";
unset output
!convert -density 300 "../img/output/oGp/iso_alka.pdf" "../img/output/oGp/iso_alka.png"
!mv "../img/output/oGp/iso_alka.pdf" "../img/output/oGp/pdf/iso_alka.pdf"

set output "../img/output/oGp/map_nly.pdf"
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
unset output
!convert -density 300 "../img/output/oGp/map_nly.pdf" "../img/output/oGp/0-map_nly.png"
!mv "../img/output/oGp/map_nly.pdf" "../img/output/oGp/pdf/map_nly.pdf"

set output "../img/output/oGp/map_isa.pdf"
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
plot "../data/output/oGp/iso_al.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "white" t "{/Symbol a}' < 0 and {/Symbol a}'' < 0", \
     0 w d lc rgb "black" t " "
unset multiplot
unset output
!convert -density 300 "../img/output/oGp/map_isa.pdf" "../img/output/oGp/1-map_isa.png"
!mv "../img/output/oGp/map_isa.pdf" "../img/output/oGp/pdf/map_isa.pdf"

set border 
set tics
set xlabel
set ylabel
set cblabel

set output "../img/output/oGp/map_isk.pdf"
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
plot "../data/output/oGp/iso_ka.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "gold" t "Re({/Symbol k}) > 0"
unset multiplot
unset output
!convert -density 300 "../img/output/oGp/map_isk.pdf" "../img/output/oGp/3-map_isk.png"
!mv "../img/output/oGp/map_isk.pdf" "../img/output/oGp/pdf/map_isk.pdf"

set border 
set tics
set xlabel
set ylabel
set cblabel

set output "../img/output/oGp/map_isak.pdf"
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
plot "../data/output/oGp/iso_al.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "white" t "{/Symbol a}' < 0 and {/Symbol a}'' < 0", \
     "../data/output/oGp/iso_ka.dat" u ($1):($2) w l dt 2 lw 2 lc rgb "gold" t "Re({/Symbol k}) > 0"
unset multiplot
unset output
!convert -density 300 "../img/output/oGp/map_isak.pdf" "../img/output/oGp/2-map_isak.png"
!mv "../img/output/oGp/map_isak.pdf" "../img/output/oGp/pdf/map_isak.pdf"
