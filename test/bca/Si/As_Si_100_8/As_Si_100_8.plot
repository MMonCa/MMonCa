set term eps enh color

set output "As_Si_100_8.eps"

set xrange[0:400]
set yrange[1e16:1e19]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 pt 6 lc rgb "black"

set label 1 "15 keV" at 30, 5e+18, 0 left norotate
set label 2 "100 keV" at 70, 2e+18, 0 left norotate

plot "as_15_8.sim" t "Test" w l ls 2, \
     "as_15_8.test" t "Ref" w l ls 1, \
     "SIMS/as15_8_30.dat" t "SIMS" w p ls 3, \
     "as_100_8.sim" notitle w l ls 2 , \
     "as_100_8.test" notitle w l ls 1, \
     "SIMS/as100_8_30.dat" notitle w p ls 3
