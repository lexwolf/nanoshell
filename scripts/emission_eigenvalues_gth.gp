set terminal pngcairo enhanced size 900,700
set output "../img/output/kap_em_multiplot.png"

set multiplot layout 2,1
set xrange [2.0:4.0]

set xlabel "ℏ{/Symbol w} (eV)"
set ylabel "ℏ{/Symbol k} (eV)"
plot "../data/output/eigenvalues_G1p1Gth.dat" u 1:10 w l t "G = 1.1xGth", \
     "../data/output/eigenvalues_G1p5Gth.dat" u 1:10 w l t "G = 1.5xGth", \
     "../data/output/eigenvalues_G1p9Gth.dat" u 1:10 w l t "G = 1.9xGth";
plot "../data/output/eigenvalues_G1p1Gth.dat" u 1:11 w l t "G = 1.1xGth", \
     "../data/output/eigenvalues_G1p5Gth.dat" u 1:11 w l t "G = 1.5xGth", \
     "../data/output/eigenvalues_G1p9Gth.dat" u 1:11 w l t "G = 1.9xGth";

unset multiplot
