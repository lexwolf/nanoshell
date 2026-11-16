# Three-gain comparison: first modulation period, omega=omega_th.
# Run from any directory by passing this script's path to gnuplot.
# Outputs: Pub_Ready/attempts/gain_switching_mechanism_composite_log.{pdf,png}
# Data columns: 1=time (ps), 6=dipole modulus, 7=population inversion.
# Gain-matched Esat values follow DATA_IMAGE_MAP.md.
# Ntarget is the common prescribed modulation, with inferred T=50 ps.
# Crossing lines delimit the first interval N > Gth/Gmax, using linear
# interpolation of actual N samples; they do not define emission onset.
# All trajectory samples within 0--50 ps are retained.

reset
script_path = ARG0
last_slash = 0
do for [i=1:strlen(script_path)] {
    if (script_path[i:i] eq '/') { last_slash = i }
}
script_dir = last_slash > 0 ? script_path[1:last_slash] : './'
cd script_dir
cd '../..'
system 'mkdir -p Pub_Ready/attempts'


array gains[3] = [1.05, 1.10, 1.15]
array tokens[3] = ['1.05Gth', '1.1Gth', '1.15Gth']
array saturation[3] = [0.325413, 0.317931, 0.310942]
# Shared logarithmic limits; initial values below 1e-15 are clipped.
intensity_min = 1e-15
intensity_max = 1e-2
array starts[3]
array ends[3]
period = 50.0
Ntarget(t) = 0.5*(1.0-cos(2.0*pi*t/period))
purple = '#8436A8'
orange = '#D97706'
blue = '#1565C0'
do for [col=1:3] {
    gain_ratio = gains[col]
    Ncrit = 1.0/gain_ratio
    data = sprintf('Data/GainSwitching/%s/omegath.dat', tokens[col])
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

    starts[col] = t_cross
    ends[col] = t_cross_end
}

set encoding utf8
set samples 2001
# Render the same layout as vector PDF and high-resolution PNG.
do for [output_kind=1:2] {
    if (output_kind == 1) {
        set terminal pdfcairo enhanced color size 7in,4.6in font 'Tinos,9' linewidth 1.2
        set output 'Pub_Ready/attempts/gain_switching_mechanism_composite_log.pdf'
    } else {
        set terminal pngcairo enhanced color size 2100,1380 font 'Tinos,27' linewidth 3.6
        set output 'Pub_Ready/attempts/gain_switching_mechanism_composite_log.png'
    }
    set border 15 linewidth 0.8
    unset grid
    set xtics in mirror scale 0.65,0.35
    set ytics in mirror scale 0.65,0.35
    set xrange [0:50]
    set xtics 10
    set mxtics 2
    set mytics 2
    set multiplot
    do for [col=1:3] {
        data = sprintf('Data/GainSwitching/%s/omegath.dat', tokens[col])
        Ncrit = 1.0/gains[col]
        Esat = saturation[col]
        set lmargin at screen (0.09+(col-1)*0.31)
        set rmargin at screen (0.355+(col-1)*0.31)
        set tmargin at screen 0.84
        set bmargin at screen 0.52
        set title sprintf('{/Tinos:Italic G}_{max} = %.2f {/Tinos:Italic G}_{th}', gains[col]) offset 0,0.2
        unset logscale y
        set mytics 2
        set yrange [0:1.08]
        set ytics 0.2
        set format y (col == 1 ? '%.1f' : '')
        set format x ''
        unset xlabel
        unset ylabel
        if (col == 1) { set ylabel 'Population inversion' offset 1,0 }
        unset key
        # One shared legend, centered above the three gain headings.
        if (col == 2) { set key at screen 0.53,0.955 center top horizontal reverse Left samplen 2 spacing 1.5 }
        set label 1 sprintf('(%s)', word('a b c',col)) at graph 0.045,0.90 front font 'Tinos:Bold'
        set arrow 1 from graph 0,first Ncrit to graph 1,first Ncrit nohead back linewidth 0.8 dashtype 3 linecolor rgb '#888888'
        set label 2 'G_{th}/G_{max}' at graph 0.045,first Ncrit offset 0,-0.9 textcolor rgb '#666666' front
        set arrow 4 from first starts[col],graph 0 to first starts[col],graph 1 nohead back linewidth 0.9 dashtype 2 linecolor rgb '#54766E'
        set arrow 5 from first ends[col],graph 0 to first ends[col],graph 1 nohead back linewidth 0.9 dashtype 2 linecolor rgb '#54766E'
        plot data using ($1<=50 ? $1 : 1/0):7 with lines linewidth 1.8 linecolor rgb purple title '{/Tinos:Italic N}(t)', \
             Ntarget(x) with lines dashtype 2 linewidth 1.6 linecolor rgb orange title '{/Tinos:Italic Ñ}(t)'
        unset title
        unset arrow 1
        unset label 2
        unset key
        set tmargin at screen 0.43
        set bmargin at screen 0.12
        set logscale y 10
        set yrange [intensity_min:intensity_max]
        set ytics (1e-14, 1e-10, 1e-6, 1e-2)
        unset mytics
        set format y '10^{%L}'
        set format x '%g'
        unset ylabel
        if (col == 1) { set ylabel '{/Tinos:Italic I}_{em} / {/Tinos:Italic I}_{sat}' offset -0.5,0 }
        set xlabel '{/Tinos:Italic t} (ps)'
        set label 1 sprintf('(%s)', word('d e f',col)) at graph 0.045,0.90 front font 'Tinos:Bold'
        plot data using ($1<=50 ? $1 : 1/0):(($6/Esat)**2) with lines linewidth 1.6 linecolor rgb blue notitle
    }
    unset multiplot
    set output
}
