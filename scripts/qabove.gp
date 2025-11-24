reset

ColorNames = "white black dark-grey red web-green web-blue dark-magenta dark-cyan dark-orange dark-yellow royalblue goldenrod dark-spring-green purple steelblue dark-red dark-chartreuse orchid aquamarine brown yellow turquoise grey0 grey10 grey20 grey30 grey40 grey50 grey60 grey70 grey grey80 grey90 grey100 light-red light-green light-blue light-magenta light-cyan light-goldenrod light-pink light-turquoise gold green dark-green spring-green forest-green sea-green blue dark-blue midnight-blue navy medium-blue skyblue cyan magenta dark-turquoise dark-pink coral light-coral orange-red salmon dark-salmon khaki dark-khaki dark-goldenrod beige olive orange violet dark-violet plum dark-plum dark-olivegreen orangered4 brown4 sienna4 orchid4 mediumpurple3 slateblue1 yellow4 sienna1 tan1 sandybrown light-salmon pink khaki1 lemonchiffon bisque honeydew slategrey seagreen antiquewhite chartreuse greenyellow gray light-gray light-grey dark-gray slategray gray0 gray10 gray20 gray30 gray40 gray50 gray60 gray70 gray80 gray90 gray100"

ColorValues = "0xffffff 0x000000 0xa0a0a0 0xff0000 0x00c000 0x0080ff 0xc000ff 0x00eeee 0xc04000 0xc8c800 0x4169e1 0xffc020 0x008040 0xc080ff 0x306080 0x8b0000 0x408000 0xff80ff 0x7fffd4 0xa52a2a 0xffff00 0x40e0d0 0x000000 0x1a1a1a 0x333333 0x4d4d4d 0x666666 0x7f7f7f 0x999999 0xb3b3b3 0xc0c0c0 0xcccccc 0xe5e5e5 0xffffff 0xf03232 0x90ee90 0xadd8e6 0xf055f0 0xe0ffff 0xeedd82 0xffb6c1 0xafeeee 0xffd700 0x00ff00 0x006400 0x00ff7f 0x228b22 0x2e8b57 0x0000ff 0x00008b 0x191970 0x000080 0x0000cd 0x87ceeb 0x00ffff 0xff00ff 0x00ced1 0xff1493 0xff7f50 0xf08080 0xff4500 0xfa8072 0xe9967a 0xf0e68c 0xbdb76b 0xb8860b 0xf5f5dc 0xa08020 0xffa500 0xee82ee 0x9400d3 0xdda0dd 0x905040 0x556b2f 0x801400 0x801414 0x804014 0x804080 0x8060c0 0x8060ff 0x808000 0xff8040 0xffa040 0xffa060 0xffa070 0xffc0c0 0xffff80 0xffffc0 0xcdb79e 0xf0fff0 0xa0b6cd 0xc1ffc1 0xcdc0b0 0x7cff40 0xa0ff20 0xbebebe 0xd3d3d3 0xd3d3d3 0xa0a0a0 0xa0b6cd 0x000000 0x1a1a1a 0x333333 0x4d4d4d 0x666666 0x7f7f7f 0x999999 0xb3b3b3 0xcccccc 0xe5e5e5 0xffffff"

myColor(c) = (idx=NaN, sum [i=1:words(ColorNames)] \
    (c eq word(ColorNames,i) ? idx=i : idx), word(ColorValues,idx))

# add transparency (alpha) a=0 to 255 or 0x00 to 0xff
myTColor(c,a) = sprintf("0x%x%s",a, myColor(c)[3:])

set term pdf color enhanced size 20cm, 20cm font ",22";
set output '../img/output/qabove.pdf';
set style line 6 lt 1 lc 0 lw 3;
set style line 7 lt 1 lc 8 lw 3;
set style line 8 dt (18,6) lc rgb myTColor("blue",180) lw 4;
set style line 9 dt (18,6) lc rgb myTColor("brown",180) lw 4;
set style line 10 lc rgb myTColor("red",80) lw 4;
set style line 11 lc rgb myTColor("black",80) lw 4;

power=10
div=10**power

set multiplot
set xrange [0:3500]
set yrange [-0.02:0.02]
set ytics -0.01, 0.01, 0.01
set key left
set size 1,0.43

set label 33 at graph 0.01, graph 0.95 "{/=20 PUMP ON}";
set object 66 rect from first 10, graph 0 to graph 1, graph 1  fc rgb "light-blue" fillstyle transparent solid;

set bmargin 4
# set lmargin 5.5
# set rmargin 7
set label 66 at graph 0.01, graph 0.1 "{/=28 (c)}"
set xlabel "Time (ps)"
set ylabel "q_2 (n. u.)"
plot "../data/output/numfunc.dat" u ($1):($6) w l ls 10 t "",\
     "../data/output/numfunc.dat" u ($1):($7) w l ls 11 t "",\
     "../data/output/anlfunc.dat" u ($1):($6) w l ls  8 t "",\
     "../data/output/anlfunc.dat" u ($1):($7) w l ls  9 t "";


unset xlabel
set format x ""
set origin 0, 0.286
set label 66 "{/=28 (b)}"
set ylabel "q_1 (n. u.)"
plot "../data/output/numfunc.dat" u ($1):($4) w l ls 10 t "",\
     "../data/output/numfunc.dat" u ($1):($5) w l ls 11 t "",\
     "../data/output/anlfunc.dat" u ($1):($4) w l ls  8 t "",\
     "../data/output/anlfunc.dat" u ($1):($5) w l ls  9 t "";
     


set origin 0, 0.572
set label 66 "{/=28 (a)}"
set ylabel "q_0 (n. u.)"
plot "../data/output/numfunc.dat" u ($1):($2) w l ls 10 t "Numerical: Real Part         ",\
     "../data/output/numfunc.dat" u ($1):($3) w l ls 11 t "Imaginary Part",\
     "../data/output/anlfunc.dat" u ($1):($2) w l ls  8 t "Analytical: Real Part         ",\
     "../data/output/anlfunc.dat" u ($1):($3) w l ls  9 t "Imaginary Part";
     
