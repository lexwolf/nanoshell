# Run from scripts/. These tails illustrate the shared detector's decisions.
reset
set terminal pngcairo size 1400,900 font ',14'
set output '../img/output/convergence_validation.png'
set label 100 'Accepted steady-state tails: paper parameters, automatic convergence' at screen 0.5,0.975 center
set multiplot layout 2,2 margins 0.09,0.98,0.08,0.90 spacing 0.12,0.19
set xlabel 'Time (ps)'
set ylabel '(I / accepted mean - 1) (%)'
set grid
set key top right font ",10" opaque
base='../data/output/convergence_validation'
do for [n=8:11] {
    set title sprintf('Raw seed E_0 = 10^{-%d}',n)
    plot sprintf('%s/automatic/%02d/accepted_tail.dat',base,n) using 1:2 with lines lc rgb 'blue' title 'Default tolerances', \
         sprintf('%s/tight/%02d/accepted_tail.dat',base,n) using 1:2 with lines dt 2 lc rgb 'dark-orange' title 'Half tolerances'
}
unset multiplot
unset output
