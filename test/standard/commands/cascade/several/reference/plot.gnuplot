set term postscript eps enh color 22
set output "defects.eps"

#set term png giant 
#set output "defects.png"

ysize=1.4
#ysize=1.
set size 0.8,ysize

set multiplot

set origin 0,0
set size 0.8,ysize/2
set yrange [10:1e5]
set xrange [70:200]
set format y "%3.0e"
set xlabel "T(K)"
set ylabel "Number of defects"
set key bottom

set logscale y
plot \
"results.txt" u 1:2  w l  t "Total", \
"results.txt" u 1:4  w lp t "I", \
"results.txt" u 1:5  w lp t "V", \
"results.txt" u 1:6  w lp t "I2", \
"results.txt" u 1:7  w lp t "V2", \
"results.txt" u 1:8  w lp t "I3", \
"results.txt" u 1:10 w lp t "I-rest", \
"results.txt" u 1:11 w lp t "V-rest


set origin 0,ysize/2
set yrange[-500:500]
unset logscale y

set xlabel "T(K)"
set ylabel "Derivative"
set key off

set arrow from 107.5,100 to 107.5,1000 nohead
set arrow from 144.5,100 to 144.5,1000 nohead
set arrow from 185.5,100 to 185.5,1000 nohead
set arrow from 278.5,100 to 278.5,1000 nohead

set arrow from 123,100 to 144,100 nohead
set arrow from 164,100 to 185,100 nohead
set arrow from 220,100 to 278,100 nohead
set arrow from 520,100 to 550,100 nohead

plot "results.txt" u ($1):(-$3) w lp t ""
