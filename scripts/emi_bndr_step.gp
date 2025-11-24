#!/bin/bash
wG=abs(ARG1);
reset
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

set key bottom
set term pdf color enhanced size 11cm, 8cm;
set lmargin screen  0.15
set rmargin screen  0.75
set tmargin screen  0.95
set tmargin screen  0.95
set bmargin screen  0.15
set output "../img/output/oGp/iso_alka.pdf"
set arrow nohead from graph 0, first wG to graph 1, first wG lc rgb "red" dt 4 lw 2 front; 
plot "../data/output/oGp/iso_al.dat" u ($1):(-$2) w filledcurves fs transparent solid 0.3 lw 2 lc rgb "black" t "{/Symbol a}' < 0 and {/Symbol a}'' < 0", \
     "../data/output/oGp/iso_ka.dat" u ($1):(-$2) w filledcurves fs transparent solid 0.3 lw 2 lc rgb "gold" t "Re({/Symbol k}) > 0";
unset output
!convert -density 300 "../img/output/oGp/iso_alka.pdf" "../img/output/oGp/iso_alka_step.png"
!mv "../img/output/oGp/iso_alka.pdf" "../img/output/oGp/pdf/iso_alka_step.pdf"
