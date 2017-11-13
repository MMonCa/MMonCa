set term postscript eps enh color 22
set output "Fe.eps"

ysize=1.4
set size 0.7,ysize

set multiplot

set origin 0,0
set size 0.7,ysize/2
set yrange [1:1e5]
set xrange [70:650]
set format y "10^{%T}"
set xlabel "T(K)"
set ylabel "Number of defects"
set key bottom

set logscale y
plot \
"results.txt" u 1:2  w l lt 1 lw 3 t "Total", \
"results.txt" u 1:4  w l lt 2 lw 3 t "I", \
"results.txt" u 1:5  w l lt 3 lw 3 t "V", \
"results.txt" u 1:6  w l lt 4 lw 3 t "I_2", \
"results.txt" u 1:7  w l lt 5 lw 3 t "V_2", \
"results.txt" u 1:8  w l lt 6 lw 3 t "I_3", \
"results.txt" u 1:9  w l lt 7 lw 3 t "Other I", \
"results.txt" u 1:10 w l lt 8 lw 3 t "Other V", \
"results.txt" u 1:11 w l lt 9 lw 3 t "<111>", \
"ChuChun-V.dat"      w p lt 3 notitle, \
"ChuChun-I.dat"      w p lt 2 notitle, \
"ChuChun-I2.dat"     w p lt 4 notitle, \
"ChuChun-V2.dat"     w p lt 5 notitle


set origin 0,ysize/2
set yrange[1:5000]
set logscale y

set xlabel "T(K)"
set ylabel "Derivative"
set key off

set arrow from 107.5,500 to 107.5,1000 head lw 2
set arrow from 144.5,500 to 144.5,1000 backhead lw 2
set arrow from 185.5,500 to 185.5,1000 backhead lw 2
set arrow from 278.5,500 to 278.5,1000 backhead lw 2
set arrow from 550.0,500 to 550.0,1000 backhead lw 2

plot "results.txt" u ($1):(-$3) w lp lw 3 t ""

unset multiplot

