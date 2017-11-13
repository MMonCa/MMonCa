set term postscript eps enh color 22
set output "test.eps"


set key at 1,1e19
set key font ",8"
set key spacing 0.5
set logscale y
# set xrange [-1:50]
# set yrange [8e19:2e20]

set xlabel "Depth (nm)"
set ylabel "Concentration (cm^-^3)" offset 0
set format y "10^{%L}"

set pointsize 1.5

set multiplot layout 2, 1

plot \
"B.init" w l lt 0 lw 3 t "B initial", \
"Bi.profile" w l lt 2 lw 3 t "Bi", \
"I.profile" w l lt 3 lw 3 t "I", \
"B_I10.init" w l lt 4 lw 3 t "B (sat) initial", \
"Bi_I10.profile" w l lt 6 lw 3 t "Bi (sat)", \
"I_I10.profile" w l lt 7 lw 3 t "I (sat)"


plot \
"B.profile" w l lt 1 lw 3 t "B", \
"B_I10.profile" w l lt 2 lw 3 t "B (sat)"

unset multiplot
