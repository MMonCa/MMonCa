set term post eps enh color
set output "nodist-test.eps"

set view map
set pm3d map
set pointsize 1.12
set palette rgbformulae 33,13,10

set size 1.1,1
set origin 0,0
set multiplot layout 2,2

set title "sxy"
splot "nodist-sxy.dat" u 2:1:3 with points palette pt 5 t ""

set title "ref sxy"
splot "data/sxy.ref" u 2:1:3 with points palette pt 5 t ""


set title "xy"
splot "nodist-xy.dat" u 2:1:3 with points palette pt 5 t ""

set title "ref"
splot "data/xy.ref" u 2:1:3 with points palette pt 5 t ""
