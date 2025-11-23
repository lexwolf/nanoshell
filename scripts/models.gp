set term pdf color enhanced size 10cm, 8cm;
set output "../img/output/models.pdf"
set xlabel "Quantity of gain G/G_{th}"
set ylabel "Emission Intensity (n.u.)"
set logscale y
plot[][] "../data/output/emission_maximum.dat" u 2:3 w l lw 2 t "|{/Symbol a}_{SS}E_0|^2" , \
         "" u 2:4 w l lw 6 dt 4 t "|p|^2 analytic", \
         "" u 2:5 w l lw 3 lc rgb "orange" t "|p|^2 numerical";
