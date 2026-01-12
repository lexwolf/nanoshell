set terminal pngcairo enhanced size 900,700
set output "../img/output/emission_eigenvalues_gth.png"

set multiplot layout 2,1
if (system("test -f ../data/output/emission_range.gp; echo $?") eq "0") { load "../data/output/emission_range.gp" } else { set xrange [2.0:4.0] }

set xlabel "ℏ{/Symbol w} (eV)"
set ylabel "ℏ{/Symbol k} (eV)"
plot "../data/output/eigenvalues_G1p1Gth.dat" u 1:10 w l t "G = 1.1xGth", \
     "../data/output/eigenvalues_G1p5Gth.dat" u 1:10 w l t "G = 1.5xGth", \
     "../data/output/eigenvalues_G1p9Gth.dat" u 1:10 w l t "G = 1.9xGth";
set key left;
plot "../data/output/eigenvalues_G1p1Gth.dat" u 1:11 w l t "G = 1.1xGth", \
     "../data/output/eigenvalues_G1p5Gth.dat" u 1:11 w l t "G = 1.5xGth", \
     "../data/output/eigenvalues_G1p9Gth.dat" u 1:11 w l t "G = 1.9xGth";

unset multiplot
