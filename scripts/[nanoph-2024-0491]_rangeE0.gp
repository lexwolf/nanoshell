# Run from scripts/, or use the companion Bash script with --plot-only.
# Based on the historical dull_rangeE0.gp, using published Fig. 1 axes.
reset
data_dir = '../data/output/[nanoph-2024-0491]_rangeE0'
parameters = data_dir.'/parameters.gp'
if (system(sprintf('test -r "%s" && echo yes', parameters)) ne 'yes') {
    print 'Error: missing completed-run parameters.gp; run the companion Bash script first.'
    exit error 'Missing normalization metadata'
}
load parameters
if (!exists('raw_seed_sweep')) { exit error 'Rerun the updated Bash script to generate the archived raw seed sweep' }
if (!(Isat > 0)) { exit error 'Isat must be positive' }

# Keep the terminal and margins from the historical dull_rangeE0.gp.
set terminal pdf color enhanced size 20cm,10cm font ',22'
set output '../img/output/[nanoph-2024-0491]_rangeE0.pdf'
set lmargin 9
set rmargin 9
set xrange [0:120]
set yrange [0:0.5]
set xtics 0,20,120
set ytics 0,0.1,0.5 nomirror
set ytics textcolor rgb 'blue'
set ylabel 'I_{em}(t)/I_{sat}' textcolor rgb 'blue'
set xlabel 'Time (ps)'
set label 44 'x10^{-3}' at graph 0,1.06 textcolor rgb 'blue'
set label 33 sprintf('{/=20 G = %.2fG_{th}}', gain_factor) at graph 0.03,0.1
set key left
set key at graph -0.03,0.97
set object 66 rect from first tpump, graph 0 to first T, graph 1 behind \
    fillcolor rgb 'light-blue' fillstyle transparent solid 0.4 border rgb 'black' linewidth 0.5

# Columns 4,5 are the complex output field; divide its squared modulus
# by Isat (0.110087 in paper mode), then display in 10^-3.
# Published legend ratios are approximate; actual raw inputs are saved per curve.
div = 1.e-3
set style line 1 linewidth 2 linecolor rgb 'blue' dashtype 1
set style line 2 linewidth 2 linecolor rgb 'blue' dashtype 2
set style line 3 linewidth 2 linecolor rgb 'blue' dashtype 3
set style line 4 linewidth 2 linecolor rgb 'blue' dashtype 4
plot for [n=8:11] sprintf('%s/numtime%02d.dat',data_dir,n) \
    using 1:(($4*$4+$5*$5)/(Isat*div)) with lines \
    linestyle (n-7) \
    title sprintf('E_0 = 10^{-%d}E_{sat}',n-1)
unset output
