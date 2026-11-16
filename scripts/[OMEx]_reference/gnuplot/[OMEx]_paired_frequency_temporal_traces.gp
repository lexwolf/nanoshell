# Typography calibrated for 5.25-inch inclusion (Cairo font scaling included): approximately 9 pt in the manuscript.
# Paired study-frequency turn-on traces at G/Gth = 1.10.
# Paths are relative to this script; rerun after replacing the input files.
# Inputs: NEW_stuff/Data/TimeBehavior/{0.9,1.11}omegath.dat.
# Intensity convention matches the existing publication figures:
# Iem/Isat = (Re(dip)^2 + Im(dip)^2)/Esat^2, Esat = 0.317931.
# Final intensity: mean of the last 100000 samples, matching TimeBehavior.ipynb.
# t_on: first crossing of 81% of each trace's own maximum intensity,
# equivalent to 90% of its maximum dipole modulus.
# linearly interpolated. This is a first-passage time, not a settling time.
# The emitted-frequency match is accepted from the supplied headers.
# All samples in the displayed 0--40 ps interval are plotted; the complete
# 200 ps records are used for the maxima and final-intensity estimates.
# Outputs: vector PDF, 300 dpi PNG, and numerical summary CSV in Pub_Ready.

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

array files[2] = ['../NEW_stuff/Data/TimeBehavior/0.9omegath.dat', '../NEW_stuff/Data/TimeBehavior/1.11omegath.dat']
array ratios[2] = [0.90, 1.11]
array colors[2] = ['#1565C0', '#D97706']
array finals[2]
array maxima[2]
array thresholds[2]
array t_on[2]
Esat = 0.317931
Nfinal = 100000
display_scale = 1e3
intensity(re,im) = (re**2+im**2)/Esat**2
set datafile commentschars '#'

do for [j=1:2] {
    stats files[j] using 1 nooutput
    count = STATS_records
    if (count < Nfinal) { print 'ERROR: too few samples for final-intensity estimate'; exit error 'Insufficient data' }
    tail_start = count-Nfinal
    stats files[j] every ::tail_start using (intensity($4,$5)) nooutput
    finals[j] = STATS_mean
    stats files[j] using (intensity($4,$5)) nooutput
    maxima[j] = STATS_max
    thresholds[j] = 0.81*maxima[j]
    target = thresholds[j]
    previous_t = 0.0
    previous_I = 0.0
    have_previous = 0
    crossing(t,y) = (hit = (!have_previous && y>=target ? t : have_previous && previous_I<target && y>=target ? previous_t+(target-previous_I)*(t-previous_t)/(y-previous_I) : 1/0), previous_t=t, previous_I=y, have_previous=1, hit)
    stats files[j] using (crossing($1,intensity($4,$5))) nooutput
    if (!exists('STATS_min')) { exit error 'No 81-percent maximum-intensity crossing found' }
    t_on[j] = STATS_min
    print sprintf('omega/omega_th=%.2f: final Iem/Isat=%.12g, t_on=%.9f ps', ratios[j], finals[j], t_on[j])
}

set print 'Pub_Ready/paired_frequency_temporal_traces_summary.csv'
print 'omega_over_omega_th,Esat,final_Iem_over_Isat,max_Iem_over_Isat,level81_Iem_over_Isat,t_on_ps,final_window_samples'
do for [j=1:2] {
    print sprintf('%.2f,%.9g,%.12g,%.12g,%.12g,%.9f,%d', ratios[j], Esat, finals[j], maxima[j], thresholds[j], t_on[j], Nfinal)
}
set print

set print $crossings
do for [j=1:2] { print sprintf('%.12g %.12g', t_on[j], display_scale*thresholds[j]) }
set print

set encoding utf8
do for [kind=1:2] {
    if (kind == 1) {
        set terminal pdfcairo enhanced color size 7in,4.5in font 'Tinos,18' linewidth 1.2
        set output 'Pub_Ready/paired_frequency_temporal_traces.pdf'
        set key font 'Tinos,18'
    } else {
        set terminal pngcairo enhanced color size 2100,1350 font 'Tinos,54' linewidth 3.6
        set output 'Pub_Ready/paired_frequency_temporal_traces.png'
        set key font 'Tinos,54'
    }
    set border 15 linewidth 0.8
    unset grid
    set xtics in mirror scale 0.65,0.35
    set ytics in mirror scale 0.65,0.35
    set xrange [0:40]
    set yrange [0:3.15]
    set xtics 5
    set ytics 0.5
    set mxtics 5
    set mytics 5
    set format y '%.1f'
    set lmargin at screen 0.16
    set rmargin at screen 0.97
    set bmargin at screen 0.19
    set tmargin at screen 0.88
    set xlabel '{/Tinos:Italic t} (ps)' offset 0,0.1
    set ylabel '{/Tinos:Italic I}_{em} / {/Tinos:Italic I}_{sat}' offset 0.5,0
    set title '{/Tinos:Italic G} = 1.10 {/Tinos:Italic G}_{th}' offset 0,0.2
    set label 1 '×10^{-3}' at graph 0,1.055 front
    set key at graph 0.97,0.08 right bottom noopaque samplen 1.5 spacing 1.2
    do for [j=1:2] {
        set arrow (j+3) from graph 0,first (display_scale*thresholds[j]) to graph 1,first (display_scale*thresholds[j]) nohead back dashtype 3 linewidth 0.8 linecolor rgb '#888888'
    }
    set label 2 '81% of peak intensity' at graph 0.97,first (display_scale*thresholds[1]) right offset 0,-1.0 textcolor rgb '#666666'
    do for [j=1:2] {
        set arrow (j+1) from first t_on[j],0 to first t_on[j],(display_scale*thresholds[j]) nohead back dashtype 3 linewidth 0.9 linecolor rgb colors[j]
    }
    plot for [j=1:2] files[j] using ($1<=40 ? $1 : 1/0):(display_scale*intensity($4,$5)) with lines linewidth 1.8 dashtype j linecolor rgb colors[j] title sprintf('ω / ω_{th} = %.2f   (t_{on} = %.2f ps)', ratios[j], t_on[j]), \
         $crossings every ::0::0 using 1:2 with points pointtype 7 pointsize 0.8 linecolor rgb colors[1] notitle, \
         $crossings every ::1::1 using 1:2 with points pointtype 7 pointsize 0.8 linecolor rgb colors[2] notitle
    set output
}
