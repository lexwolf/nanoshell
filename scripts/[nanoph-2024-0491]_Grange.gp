reset
repo = (system("test -d scripts -a -d data -a -d img; echo $?") eq "0") ? "." : ".."
#  omin = 2.; omax = 3.4;
omin = 2.6; omax = 3.;

# SETTING THE VISIBLE SPECTRUM IMAGE AT THE BOTTOM
set samples 200
set isosamples 2
k=omax-omin
set cbrange [omin:omax]
r(x)=x<2.1?1:x<2.4?-(x-2.4)/(2.4-2.1):x<2.8?0:x<3.3?(x-2.8)/(3.3-2.8):1
g(x)=x<1.9?0:x<2.1?(x-1.9)/(2.1-1.9):x<2.5?1:x<2.8?-(x-2.8)/(2.8-2.5):0
b(x)=x>2.5?1:x>2.4?-(x-2.4)/(2.4-2.5):0
f(x)=x<1.5?0:x<1.8?0.3+0.7*(x-1.6)/(1.8-1.6):x<3.?1:x<3.41?0.3+0.7*(3.3-x)/(3.3-3.):0
set palette functions f(k*gray+omin)*r(k*gray+omin),g(k*gray+omin),f(k*gray+omin)*b(k*gray+omin)
set palette functions f(k*gray+omin)*r(k*gray+omin),g(k*gray+omin),f(k*gray+omin)*b(k*gray+omin)
unset colorbox
# DONE

# DEFINING COLORS #
color125 = "dark-salmon"
color150 = "purple"
color175 = "skyblue"
# DONE

omeB = 2.8121972

set size .95,.95
set term pdf enhanced size 23cm, 20cm font ",26"
set output repo."/img/output/[nanoph-2024-0491]_Grange.pdf"
set multiplot

set lmargin at screen 0.16
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.48
set ytics nomirror tc rgb "red"
set y2tics tc rgb "blue"
set key left
set xlabel "G/G_{th}"
set yrange [0:0.1]

power=-2
div=10**power
set label 132 sprintf("x10^{%d}", power) at graph 0.95, 1.07 tc rgb "blue"
set label 55 at graph 0.02, graph 0.1 "{/=32 (d)}"

set ylabel "{/Symbol D}h{/Symbol w}_{em} (eV)" tc rgb "red"
set y2label "I_{max}/I_{sat}" tc rgb "blue"

set arrow 10 nohead from first 1.25, graph 0. to first 1.25, graph 1 lw 4 dt 2
set label 77 at first 1.14, graph 0.94 "(a)"
set arrow 20 nohead from first 1.5, graph 0. to first 1.5, graph 1 lw 4 dt 2
set label 78 at first 1.39, graph 0.94 "(b)"
set arrow 30 nohead from first 1.75, graph 0. to first 1.75, graph 1 lw 4 dt 2
set label 79 at first 1.65, graph 0.94 "(c)"

set y2tics 0,1.5
plot repo."/data/output/Grange/emission_maximum.dat" u ($2):($5/($8*div)) axis x1y2 w filledcurves y2 fs transparent solid 0.50 lw 1 lc palette cb omeB t "I_{max}/I_{sat}", \
     "" u ($2):($5/($8*div)) axis x1y2 w l lw 4 lc rgb "blue" t "", \
     "" u ($2):($6) w l lw 4 lc rgb "red" t "{/Symbol D}h{/Symbol w}_{em}"
unset label 77
unset label 78
unset label 79

set ytics mirror tc rgb "blue"
set ytics 0,1.5
unset arrow
unset object
unset y2tics

set label 132 at graph 0., 1.07 tc rgb "blue"
set ylabel "I_{em}/I_{sat}" tc rgb "blue"

unset y2label
set xlabel "h{/Symbol w}_{em} (eV)"
set label 66 "G = 1.25 G_{th}" at graph 0.05, graph 0.9
Isat = 0.0889504

set lmargin at screen 0.16
set rmargin at screen 0.39
set bmargin at screen 0.63
set tmargin at screen 0.95

set xtics 2.6,0.1,2.9
set label 55 "{/=32 (a)}"
plot [omin:omax][0:4.5] repo."/data/output/Grange/oG1.25p/ome_p3.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB t ""

set lmargin at screen 0.39
set rmargin at screen 0.63
set bmargin at screen 0.63
set tmargin at screen 0.95

set label 66 "G = 1.5 G_{th}" at graph 0.05, graph 0.9
Isat = 0.0741253
unset ylabel
set format y ""
unset label 132
set label 55 "{/=32 (b)}"
set xtics 2.6,0.1,2.9
plot [omin:omax][0:4.5] repo."/data/output/Grange/oG1.5p/ome_p3.dat" u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB t ""
set y2label "I_{em}/I_{sat}" tc rgb "blue"

set lmargin at screen 0.63
set rmargin at screen 0.85
set bmargin at screen 0.63
set tmargin at screen 0.95

set label 66 "G = 1.75 G_{th}" at graph 0.05, graph 0.9
Isat = 0.063536
set label 55 "{/=32 (c)}"
set label 132 sprintf("x10^{%d}", power) at graph 0.82, 1.07 tc rgb "blue"
set y2tics auto tc rgb "blue"
set y2tics 0,1.5
set xtics 2.6,0.1,3.
set y2range [0:4.5]
plot [omin:omax][0:4.5] repo."/data/output/Grange/oG1.75p/ome_p3.dat" u ($1):($2/(div*Isat)) axis x1y2 w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omeB t ""

unset multiplot
unset output
