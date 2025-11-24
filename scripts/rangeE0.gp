reset
# DEFINING COMMAND TO READ FROM files
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
# DONE

# READING TOTAL TIME AND TIME_PUMP_ON (ome T tpump) FROM "../data/input/time.dat"
file="../data/input/time.dat" ;
row=1 ; col=1
T = abs(at(file,row,col))
row=1 ; col=2
tpump = at(file,row,col)
# DONE

ColorNames = "white black dark-grey red web-green web-blue dark-dark-violet dark-cyan dark-orange dark-yellow royalblue goldenrod dark-spring-green purple steelblue dark-red dark-chartreuse orchid aquamarine brown yellow turquoise grey0 grey10 grey20 grey30 grey40 grey50 grey60 grey70 grey grey80 grey90 grey100 light-red light-green light-blue light-dark-violet light-cyan light-goldenrod light-pink light-turquoise gold green dark-green spring-green forest-green sea-green blue dark-blue midnight-blue navy medium-blue skyblue cyan dark-violet dark-turquoise dark-pink coral light-coral orange-red salmon dark-salmon khaki dark-khaki dark-goldenrod beige olive orange violet dark-violet plum dark-plum dark-olivegreen orangered4 brown4 sienna4 orchid4 mediumpurple3 slateblue1 yellow4 sienna1 tan1 sandybrown light-salmon pink khaki1 lemonchiffon bisque honeydew slategrey seagreen antiquewhite chartreuse greenyellow gray light-gray light-grey dark-gray slategray gray0 gray10 gray20 gray30 gray40 gray50 gray60 gray70 gray80 gray90 gray100"

ColorValues = "0xffffff 0x000000 0xa0a0a0 0xff0000 0x00c000 0x0080ff 0xc000ff 0x00eeee 0xc04000 0xc8c800 0x4169e1 0xffc020 0x008040 0xc080ff 0x306080 0x8b0000 0x408000 0xff80ff 0x7fffd4 0xa52a2a 0xffff00 0x40e0d0 0x000000 0x1a1a1a 0x333333 0x4d4d4d 0x666666 0x7f7f7f 0x999999 0xb3b3b3 0xc0c0c0 0xcccccc 0xe5e5e5 0xffffff 0xf03232 0x90ee90 0xadd8e6 0xf055f0 0xe0ffff 0xeedd82 0xffb6c1 0xafeeee 0xffd700 0x00ff00 0x006400 0x00ff7f 0x228b22 0x2e8b57 0x0000ff 0x00008b 0x191970 0x000080 0x0000cd 0x87ceeb 0x00ffff 0xff00ff 0x00ced1 0xff1493 0xff7f50 0xf08080 0xff4500 0xfa8072 0xe9967a 0xf0e68c 0xbdb76b 0xb8860b 0xf5f5dc 0xa08020 0xffa500 0xee82ee 0x9400d3 0xdda0dd 0x905040 0x556b2f 0x801400 0x801414 0x804014 0x804080 0x8060c0 0x8060ff 0x808000 0xff8040 0xffa040 0xffa060 0xffa070 0xffc0c0 0xffff80 0xffffc0 0xcdb79e 0xf0fff0 0xa0b6cd 0xc1ffc1 0xcdc0b0 0x7cff40 0xa0ff20 0xbebebe 0xd3d3d3 0xd3d3d3 0xa0a0a0 0xa0b6cd 0x000000 0x1a1a1a 0x333333 0x4d4d4d 0x666666 0x7f7f7f 0x999999 0xb3b3b3 0xcccccc 0xe5e5e5 0xffffff"

 omin = 2.; omax = 3.4;
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

omeB = 2.8121972
 
set term pdf color enhanced size 20cm, 26cm font ",22";
set output '../img/output/rangeE0.pdf';
set style line 6 lt 1 lc 0 lw 3;
set style line 7 lt 1 lc 8 lw 3;
set style line 8 lc rgb "purple" lw 4;
set style line 9 lc rgb "black" lw 4;
set style line 10 lt 1.5 lc rgb "dark-violet" lw 3;
set style line 11 lt 1.5 lc rgb "dark-green" lw 2;
set style line 12 lt 1.5 lc rgb "dark-violet" lw 4;

set style fill transparent solid 0.4

power=-3
div=10**power

# set xrange [0:800]
set key left


set label 33 at graph 0.03, graph 0.1 "{/=20 G = 1.01G_{th}}";
set object 66 rect from first tpump, graph 0 to first T, graph 1  fc rgb "light-blue" fillstyle transparent solid;

set lmargin 9
set rmargin 9

set ytics nomirror 0,0.01
set ytics tc rgb "blue"
set ylabel "Emission Intensity" tc rgb "blue" offset 1,-12
set y2tics tc rgb "dark-violet"
set y2label "Population Inversion" tc rgb "dark-violet" offset 0,-12
set label 44 sprintf("x10^{%d}", power-2) at graph 0, 1.05 tc "blue"

set y2tics 0.98, 0.01
set y2range [.98:1]
set yrange [0:0.05]
set ytics 0.01, 0.01
set format y2 "%.2f"
set format x ""
ts=0.03
sp=0.
mx=1-sp-ts
ht=0.22
dh=sp+ht
set multiplot
set label 22 at graph 0.024, graph 0.9 font "/,32" "(a)"
set label 23 at graph 0.03, graph 0.75 "E_0 = 10^{-9}E_{sat}"
set tmargin at screen mx; set bmargin at screen mx-ht
set y2tics 0.99, 0.01
# plot sin(x) lc rgb "blue"
plot "../data/output/numtime08.dat" u ($1):(($4*$4+$5*$5)/div) w filledcurve x1 fs transparent solid 0.50 lw 1 lc palette cb omeB t "",\
     "../data/output/numtime08.dat" u ($1):(($4*$4+$5*$5)/div) w l lt 1 lw 2 lc rgb "blue" t "",\
     "../data/output/numtime08.dat" u ($1):($6) axis x1y2 w l ls 12 t "";
unset ylabel
unset y2label
unset label 44
set label 22 "(b)"
set label 23 "E_0 = 10^{-10}E_{sat}"
set tmargin at screen mx-dh; set bmargin at screen mx-dh-ht
# plot sin(x) lc rgb "blue"
plot "../data/output/numtime09.dat" u ($1):(($4*$4+$5*$5)/div) w filledcurve x1 fs transparent solid 0.50 lw 1 lc palette cb omeB t "",\
     "../data/output/numtime09.dat" u ($1):(($4*$4+$5*$5)/div) w l lt 1 lw 2 lc rgb "blue" t "",\
     "../data/output/numtime09.dat" u ($1):($6) axis x1y2 w l ls 12 t "";
set label 22 "(c)"
set label 23 "E_0 = 10^{-11}E_{sat}"
set tmargin at screen mx-2*dh; set bmargin at screen mx-2*dh-ht
# plot sin(x) lc rgb "blue"
plot "../data/output/numtime10.dat" u ($1):(($4*$4+$5*$5)/div) w filledcurve x1 fs transparent solid 0.50 lw 1 lc palette cb omeB t "",\
     "../data/output/numtime10.dat" u ($1):(($4*$4+$5*$5)/div) w l lt 1 lw 2 lc rgb "blue" t "",\
     "../data/output/numtime10.dat" u ($1):($6) axis x1y2 w l ls 12 t "";
set label 22 "(d)"
set label 23 "E_0 = 10^{-12}E_{sat}"
set tmargin at screen mx-3*dh; set bmargin at screen mx-3*dh-ht
set xlabel "Time (ps)"
set format x "%g"
set y2tics 0.98, 0.01
set ytics 0., 0.01
set key at graph -0.04, graph 0.42
# plot sin(x) lc rgb "blue"
plot "../data/output/numtime11.dat" u ($1):(($4*$4+$5*$5)/div) w filledcurve x1 fs transparent solid 0.50 lw 1 lc palette cb omeB t "I_{em}(t)/I_{sat}",\
     "../data/output/numtime11.dat" u ($1):(($4*$4+$5*$5)/div) w l lt 1 lw 2 lc rgb "blue" t "",\
     "../data/output/numtime11.dat" u ($1):($6) axis x1y2 w l ls 12 t "N(t)";
unset multiplot
unset output
