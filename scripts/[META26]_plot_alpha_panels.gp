set terminal pdfcairo enhanced color size 34cm,8cm font "Arial,18"

repo = (system("test -d scripts -a -d data -a -d img; echo $?") eq "0") ? "." : ".."
set output repo."/img/output/alpha_panels.pdf"

data_dir = repo."/data/output/alpha_panels"
meta = data_dir."/metadata.dat"

f0 = data_dir."/stationary_G_0.00Gth.dat"
f1 = data_dir."/stationary_G_0.40Gth.dat"
f2 = data_dir."/stationary_G_0.90Gth.dat"
f3 = data_dir."/stationary_G_1.00Gth.dat"
f4 = data_dir."/stationary_G_1.04Gth.dat"
f5 = data_dir."/stationary_G_1.60Gth.dat"

ome_sp = real(system(sprintf("awk -F= '/^ome_sp/ {gsub(/[[:space:]]/, \"\", $2); print $2}' %s", meta)))

xmin = 1e300
xmax = -1e300

stats f0 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)
stats f1 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)
stats f2 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)
stats f3 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)
stats f4 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)
stats f5 using 1 nooutput; xmin = (STATS_min < xmin ? STATS_min : xmin); xmax = (STATS_max > xmax ? STATS_max : xmax)

half_max_cmd = "awk 'NF >= 3 && $1 !~ /^#/ { n++; x[n]=$1; y[n]=$3; if (n == 1 || y[n] > ymax) { ymax=y[n]; imax=n } } END { if (n < 3 || imax <= 1 || imax >= n) exit 1; half=0.5*ymax; found_left=0; found_right=0; for (i=1; i<imax; i++) { if ((y[i]-half)*(y[i+1]-half) <= 0 && y[i] != y[i+1]) { omega_left=x[i]+(half-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); found_left=1 } } for (i=imax; i<n; i++) { if ((y[i]-half)*(y[i+1]-half) <= 0 && y[i] != y[i+1]) { omega_right=x[i]+(half-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); found_right=1; break } } if (!found_left || !found_right) exit 1; delta=(x[imax]-omega_left > omega_right-x[imax] ? x[imax]-omega_left : omega_right-x[imax]); if (delta <= 0) exit 1; printf \"%.12g %.12g\\n\", x[imax]-8*delta, x[imax]+8*delta }' ".f0
half_max_range = system(half_max_cmd)
if (words(half_max_range) >= 2) {
    xmin = real(word(half_max_range, 1))
    xmax = real(word(half_max_range, 2))
}

stats f1 using (abs($2)) nooutput
alpha_abs_real_MAX = STATS_max
stats f1 using (abs($3)) nooutput
alpha_abs_imag_MAX = STATS_max
alpha_abs_MAX_raw = (alpha_abs_real_MAX > alpha_abs_imag_MAX ? alpha_abs_real_MAX : alpha_abs_imag_MAX)
alpha_abs_pow = 10.0 ** floor(log10(alpha_abs_MAX_raw))
alpha_abs_MAX = ceil(alpha_abs_MAX_raw / alpha_abs_pow) * alpha_abs_pow
if (alpha_abs_MAX == 0) alpha_abs_MAX = 1

set xrange [xmin:xmax]
set yrange [-1.05 * alpha_abs_MAX:1.05 * alpha_abs_MAX]

set style line 1 lc rgb "black" lw 5
set style line 2 lc rgb "#e9551a" lw 5
set style line 3 lc rgb "#666666" lw 1.5 dt 2
title_re = "{/Symbol a}'/4{/Symbol p}a^3"
title_im = "{/Symbol a}''/4{/Symbol p}a^3"

set border lw 1.8
set tics out nomirror
set xtics 0.2
set mxtics 2
set format y "%.0f"
set xlabel "h{/Symbol w} (eV)"
set ylabel "Polarizability ({/Symbol a}/4{/Symbol p}a^3)"
set grid noxtics ytics lc rgb "#dddddd" lw 0.5
set bmargin 4
set multiplot layout 1,6 margins 0.075,0.995,0.30,0.95 spacing 0.018,0

unset key
set arrow 1 from ome_sp, graph 0 to ome_sp, graph 1 nohead ls 3 front
set label 2 "G = 0 Gth" at graph 0.01,0.08 front font ",24"
plot f0 using 1:2 with lines ls 1 title title_re, \
     f0 using 1:3 with lines ls 2 title title_im

unset ylabel
unset ytics
set label 2 "G = 0.4 Gth" at graph 0.01,0.08 front font ",24"
plot f1 using 1:2 with lines ls 1 title title_re, \
     f1 using 1:3 with lines ls 2 title title_im

set label 2 "G = 0.9 Gth" at graph 0.01,0.08 front font ",24"
plot f2 using 1:2 with lines ls 1 title title_re, \
     f2 using 1:3 with lines ls 2 title title_im

set label 2 "G = 1.0 Gth" at graph 0.01,0.08 front font ",24"
plot f3 using 1:2 with lines ls 1 title title_re, \
     f3 using 1:3 with lines ls 2 title title_im

set label 2 "G = 1.04 Gth" at graph 0.01,0.08 front font ",24"
plot f4 using 1:2 with lines ls 1 title title_re, \
     f4 using 1:3 with lines ls 2 title title_im

set key right top spacing 0.85 samplen 1.8 font ",24"
set label 2 "G = 1.6 Gth" at graph 0.01,0.08 front font ",24"
plot f5 using 1:2 with lines ls 1 title title_re, \
     f5 using 1:3 with lines ls 2 title title_im

unset multiplot
