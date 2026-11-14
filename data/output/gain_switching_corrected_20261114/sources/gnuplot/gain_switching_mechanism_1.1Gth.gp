# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Two-panel gain-switching mechanism at Gmax/Gth=1.1, omega=omega_th.
# Run: gnuplot gnuplot/gain_switching_mechanism_1.1Gth.gp (from Gain_Switching).
# Data: ordinary threshold trajectory, columns 1=t(ps), 6=|dip|, 7=N.
# Target inversion: Data_Description.pdf Fig.13, (1-cos(2*pi*t/T))/2.
# T=50 ps is inferred from the repeated trajectory (pulse spacing ~50 ps),
# consistent with Fig.24; the simulation's parameter file is not supplied.
# Ncrit=Gth/Gmax is the static threshold reference, not a prediction of the
# exact transient emission onset. Esat is the gain-matched value .317931.
# All samples in the displayed interval are retained (no downsampling).

reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '..'
system 'mkdir -p Pub_Ready'

gain_ratio = 1.1
Esat = 0.317931
period = 50.0
Ncrit = 1.0/gain_ratio
Ntarget(t) = 0.5*(1.0-cos(2.0*pi*t/period))
data = 'Data/GainSwitching/1.1Gth/omegath.dat'
purple = '#8436A8'
orange = '#D97706'
blue = '#1565C0'

# Interpolate the first upward and subsequent downward crossings of Ncrit.
# These delimit the first interval above the static threshold.
set datafile commentschars '#'
previous_t = 0.0
previous_N = 0.0
crossing(t,n) = (hit = (previous_N<Ncrit && n>=Ncrit ? previous_t+(Ncrit-previous_N)*(t-previous_t)/(n-previous_N) : 1/0), previous_t=t, previous_N=n, hit)
stats data using (crossing($1,$7)) nooutput
t_cross = STATS_min
previous_t = 0.0
previous_N = 0.0
crossing_down(t,n) = (hit = (t>t_cross && previous_N>=Ncrit && n<Ncrit ? previous_t+(Ncrit-previous_N)*(t-previous_t)/(n-previous_N) : 1/0), previous_t=t, previous_N=n, hit)
stats data using (crossing_down($1,$7)) nooutput
t_cross_end = STATS_min
print sprintf('Above-threshold interval: %.6f to %.6f ps (duration %.6f ps)', t_cross, t_cross_end, t_cross_end-t_cross)

set encoding utf8
set terminal pdfcairo enhanced color size 7in,6in font 'Tinos,18' linewidth 1.5
set output 'Pub_Ready/gain_switching_mechanism_1.1Gth.pdf'
set datafile commentschars '#'
set samples 2001
set border 15 linewidth 1
unset grid
set xtics in mirror scale 0.65,0.35
set ytics in mirror scale 0.65,0.35
unset x2tics
unset y2tics
set xrange [0:100]
set xtics 20
set mxtics 4
set mytics 2
set lmargin at screen 0.17
set rmargin at screen 0.975

set multiplot
set tmargin at screen 0.92
set bmargin at screen 0.54
set title '{/Tinos:Italic G}_{max} = 1.1 {/Tinos:Italic G}_{th},  ω = ω_{th}' font 'Tinos,18' offset 0,-0.4
set yrange [0:1.08]
set ytics 0.2
set format y '%.1f'
set format x ''
set ylabel 'Population inversion' font 'Tinos,18' offset 0.5,0
set key right top vertical reverse Left noopaque samplen 2 spacing 1.5
set label 1 '(a)' at graph 0.025,0.92 front font 'Tinos:Bold,18'
# Use a precomposed Unicode tilde-N for reliable Cairo legend rendering.
set arrow 1 from graph 0,first Ncrit to graph 1,first Ncrit nohead back linewidth 0.8 dashtype 3 linecolor rgb '#888888'
set label 2 'G_{th}/G_{max}' at graph 0.025,first Ncrit offset 0,-0.9 textcolor rgb '#666666' font 'Tinos,18' front
set arrow 4 from first t_cross,graph 0 to first t_cross,graph 1 nohead back linewidth 0.9 dashtype 2 linecolor rgb '#54766E'
set arrow 5 from first t_cross_end,graph 0 to first t_cross_end,graph 1 nohead back linewidth 0.9 dashtype 2 linecolor rgb '#54766E'
plot data using ($1<=100 ? $1 : 1/0):7 with lines linewidth 1.8 linecolor rgb purple title '{/Tinos:Italic N}(t)', \
     Ntarget(x) with lines dashtype 2 linewidth 1.6 linecolor rgb orange title '{/Tinos:Italic Ñ}(t)'

unset title
unset arrow 1
unset label 2
unset key
set tmargin at screen 0.46
set bmargin at screen 0.15
set yrange [0:1.3]
set ytics 0.4
set format y '%.1f'
set label 6 '×10^{-3}' at graph 0,1.05 front font 'Tinos,18'
set format x '%g'
set ylabel '{/Tinos:Italic I}_{em} / {/Tinos:Italic I}_{sat}' font 'Tinos,18' offset 0.5,0
set xlabel '{/Tinos:Italic t} (ps)' font 'Tinos,18' offset 0,0.1
set label 1 '(b)' at graph 0.025,0.92 front font 'Tinos:Bold,18'
plot data using ($1<=100 ? $1 : 1/0):(1e3*($6/Esat)**2) with filledcurves y=0 fillstyle transparent solid 0.3 noborder linecolor rgb blue notitle, \
     data using ($1<=100 ? $1 : 1/0):(1e3*($6/Esat)**2) with lines linewidth 1.6 linecolor rgb blue notitle

unset multiplot
set output
