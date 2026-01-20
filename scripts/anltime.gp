#!/bin/bash
reset 
# DEFINING COMMAND TO READ FROM files
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
# DONE

# READING Gth FROM "../data/output/frohlich.dat"
file="../data/output/frohlich.dat"
row=1 ; col=2
Gth = at(file,row,col)
# DONE

# CALCULATING GAIN RANGE
rGth = 1.1 * Gth
# DONE

# READING SPECTRUM RANGE (omin omax),  ome_g AND G FROM "../data/input/nanosphere_eV.dat"
file="../data/input/nanosphere_eV.dat" ;
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)
row=1 ; col=3
ome_g = at(file,row,col)
row=1 ; col=4
G = abs(at(file,row,col))
row=1 ; col=12
rho = abs(at(file,row,col))
# DONE

# Define functions to calculate color values
r(x) = x < 2.1 ? 1 : x < 2.4 ? -(x-2.4)/(2.4-2.1) : x < 2.8 ? 0 : x < 3.3 ? (x-2.8)/(3.3-2.8) : 1
g(x) = x < 1.9 ? 0 : x < 2.1 ? (x-1.9)/(2.1-1.9) : x < 2.5 ? 1 : x < 2.8 ? -(x-2.8)/(2.8-2.5) : 0
b(x) = x > 2.5 ? 1 : x > 2.4 ? -(x-2.4)/(2.4-2.5) : 0
f(x) = x < 1.5 ? 0 : x < 1.8 ? 0.3 + 0.7*(x-1.8)/(1.8-1.6) : x < 3. ? 1 : x < 3.41 ? 0.3 + 0.7*(3.3-x)/(3.3-3.) : 0

# Define the color values
r_value = f(ome_g) * r(ome_g)
g_value = g(ome_g)
b_value = f(ome_g) * b(ome_g)

# Format the RGB values to hexadecimal
rgb_color = sprintf("#%02x%02x%02x", int(255 * r_value), int(255 * g_value), int(255 * b_value))

# Define the new palette
set palette defined (0 "white", 1 rgb_color)

set cbrange [0:1.1]

# Function to set the background for each plot
set_background(file, ymax) = sprintf("unset title; unset xlabel; set pm3d map; unset xtics; unset ytics; unset border; splot '%s' u ($1):(%f*$2):($3 / Gth) t ''; unset pm3d; set border; set ytics; set xtics;", file, ymax)

# Function to set margins
set_margins(lmar, rmar, tmar, bmar) = sprintf("set lmargin at screen %.4f; set rmargin at screen %.4f; set tmargin at screen %.4f; set bmargin at screen %.4f;", lmar, rmar, tmar, bmar)

# READING FREQUENCY (ome) FROM "../data/output/omega.dat"
file="../data/output/omega.dat" ;
row=1 ; col=1
omega = abs(at(file,row,col))
# DONE

# READING TOTAL TIME AND TIME_PUMP_ON (ome T tpump) FROM "../data/input/time.dat"
file="../data/input/time.dat" ;
row=1 ; col=1
T = abs(at(file,row,col))
row=1 ; col=2
tpump = at(file,row,col)
# DONE

# READING STATIONARY QUASI STATIC ALPHA (Re(alph_QS), Im(alph_QS)) AND Im(eps2) FROM "../data/output/alpha.dat"
file="../data/output/alpha.dat" ;
row=1 ; col=1
reaQS = at(file,row,col)
row=1 ; col=2
imaQS = at(file,row,col)
row=1 ; col=6
ieps1 = abs(at(file,row,col))
# DONE

# Define margins
left_margin = 0.08
right_margin = 0.9
top_margin = 0.88
bottom_margin = 0.2

# Calculate plot width and height
plot_width = (right_margin - left_margin) / 3
plot_height = (top_margin - bottom_margin)

stats "../data/output/stationary.dat" using 2 name "Y2" nooutput
stats "../data/output/stationary.dat" using 3 name "Y3" nooutput

# Calcola il valore massimo tra le due colonne
ymax = (Y2_max > Y3_max) ? Y2_max : Y3_max
ymin = (Y2_min < Y3_min) ? Y2_min : Y3_min

# Imposta il terminale e la gamma y
set term unknown
if (ymax > 200) {
    ymin = -200
    ymax = 200
    } else {
    ymin = 1.1*ymin
    ymax = 1.1*ymax
    }

set yrange [ymin:ymax]
set y2range [ymin:ymax]
    
set xrange [0:T];
set term pdf color enhanced size 20cm, 8cm font ",16";
set output '../img/output/intime.pdf';
set style line 6 lt 1 lc 0 lw 4;
set style line 7 lt 1 lc 8 lw 4;
set size 2,1.;

Gpal = ieps1 / rGth;
Gnor = G / Gth;

set title sprintf("G = %.1f × G_{th} —  {/Symbol r} = %.1f", Gnor, rho);
set cblabel "|{/Symbol e}@{''}_2({/Symbol w})| / G_{th}"

set multiplot layout 2, 2 margins screen left_margin, screen right_margin, screen bottom_margin, screen top_margin
eval set_margins(left_margin, left_margin + 1.79 * plot_width, top_margin, bottom_margin)
set ylabel "Normalized Polarizability {/Symbol a}/(4{/Symbol p}a^3)"

# set grid ytics
set label 22 at graph 0.02, 0.08 "{/=24 (a)}";

if (G!=0) { 
    set label 33 at first 0.22, graph 0.95 "{/=14 PUMP ON}";
    set object 66 rect from first tpump, graph 0 to first T, graph 1  fc palette frac Gpal fs solid back;
}
set xlabel "Time (ps)"

plot "../data/output/anltime.dat" u 1:2 w l ls 6 t "Re[{/Symbol a}]/(4{/Symbol p}a^3)",\
     "../data/output/anltime.dat" u 1:3 w l ls 7 t "Im[{/Symbol a}]/(4{/Symbol p}a^3)";
     
dome=(omax-omin)/8.
dscale=10**floor(log10(dome))
dbase=dome/dscale
dome=(dbase < 2.5 ? dscale : (dbase <= 7.5 ? 5*dscale : 10*dscale))

set xtics omin+dome, dome
set xrange [omin:omax];
set label 22 at graph 0.02, 0.08 "{/=24 (b)}";
unset label 33
unset ylabel
set format y ""

set origin 0.52,0
set rmargin 17


eval set_margins(left_margin + 1.81 * plot_width, right_margin, top_margin, bottom_margin)
eval set_background("../data/output/background.dat", 1 * ymax)
eval set_margins(left_margin + 1.81 * plot_width, right_margin, top_margin, bottom_margin)
set arrow 22 nohead from first omega, graph 0 to first omega, graph 1 lc rgb "dark-green" lw 2 front
set arrow 23 nohead from graph 0, first reaQS to first omega, first reaQS lc 0 lw 2 dt 2 front
set arrow 24 nohead from graph 0, first imaQS to first omega, first imaQS ls 7 lw 2 dt 2 front
set title sprintf("ℏ{/Symbol w} = %.2f",omega) tc rgb "dark-green";
set xlabel "Energy (eV)"

plot "../data/output/stationary.dat" u 1:2 w l ls 6 t "",\
     "../data/output/stationary.dat" u 1:3 w l ls 7 t ""; 

unset multiplot
unset term
! convert -density 200 ../img/output/intime.pdf ../img/output/intime.png
