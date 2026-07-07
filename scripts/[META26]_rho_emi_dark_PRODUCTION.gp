reset
omin=2.8
omax=4.6
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
r_list="0.8 0.7 0.6 0.5 0.4"
n_r=words(r_list)
xmax_list="4.0616 4.00718 3.95415 3.90505 3.86317"
ymax_list="0.000111855 9.62592e-05 6.56003e-05 3.49817e-05 1.35141e-05"
dx_list="0.044 .132 .220 .308 .396"

# Editable production layout, in plot coordinates.
# Order follows r_list. Change these values, then rerun:
#   gnuplot [nanoph-2024-0491]_rho_emi.gp
sketch_x_list="3.85 3.79 3.75 3.7 3.66"
sketch_y_list="0.01086478366 0.0085 0.0058 0.0025 0.0002"
label_x_list="3.6 3.54 3.5 3.45 3.41"
label_y_list="0.01195126366 0.0095 0.0068 0.0035 0.0012"

omegaB(i) = at(sprintf("../data/output/rho/omeB-%s.dat", word(r_list,i)),1,1)
xmax(i) = real(word(xmax_list,i))
ymax_raw(i) = real(word(ymax_list,i))
dx(i) = real(word(dx_list,i))
sketch_x(i) = real(word(sketch_x_list,i))
sketch_y(i) = real(word(sketch_y_list,i))
label_x(i) = real(word(label_x_list,i))
label_y(i) = real(word(label_y_list,i))

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

set term pdf color enhanced size 10cm, 8cm;
set output "../img/output/[META26]_rho_dark_PRODUCTION.pdf"

set multiplot
# PLOTTING THE VISIBLE SPECTRUM
set origin 0,0.03
set size 0.991,0.25
set pm3d map
unset ytics
set lmargin at screen 0.16
set rmargin at screen 0.97
set bmargin at screen 0.15
set tmargin at screen 0.2
set xlabel "ℏ{/Symbol w}_{em} (eV)" offset 0,-0.3
set xtics offset 0,-0.5
splot[omin:omax] x t ""
# DONE

Isat=0.0926567

unset xlabel

set bmargin at screen 0.2
set tmargin at screen 0.94
unset xtics

set yrange [:0.2]
set ytics 0, 0.004

power=-1
div=10**power

scaleY(y)=y/(div*Isat)

# Determine dynamic y-range and place thumbnails/labels near peaks
max_plot_y = 0
do for [i=1:n_r] {
    py = scaleY(ymax_raw(i))
    if (py > max_plot_y) { max_plot_y = py }
}
if (max_plot_y <= 0) { max_plot_y = 1 }
set yrange [0:max_plot_y*1.1]

set for [i=1:n_r] pixmap (3+i) sprintf("../data/output/rho/%s.png", word(r_list,i)) at first sketch_x(i), first sketch_y(i) width screen 0.08
set for [i=1:n_r] label sprintf("{/Symbol r} = %s", word(r_list,i)) at first label_x(i), first label_y(i) left

set ylabel "I_{em}/I_{sat}"
set xrange [omin:omax]
set label sprintf("x10^{%d}", power) at graph 0, 1.03

plot for [i=n_r:1:-1] sprintf("../data/output/rho/%s.dat", word(r_list,i)) u ($1):($2/(div*Isat)) w filledcurves fs transparent solid 0.50 lw 1 lc palette cb omegaB(i) t ""
unset multiplot
unset output
