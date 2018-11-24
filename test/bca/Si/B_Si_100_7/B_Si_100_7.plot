set term eps enh color

set output "B_Si_100_7.eps"

set xrange[0:4000]
set yrange[1e16:2e18]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

set label 1 "15 keV" at 200, 1.5e+18, 0 left norotate
set label 2 "80 keV" at 400, 1e+18, 0 left norotate
set label 3 "280 keV" at 800, 5e+17, 0 left norotate
set label 4 "700 keV" at 1600, 4e+17, 0 left norotate
set label 5 "2400 keV" at 3000, 5e+17, 0 left norotate

plot "b_15_7.sim" t "Test" w l ls 2, \
     "b_15_7.test" t "Ref" w l ls 1, \
     "SIMS/b15_7.dat" t "SIMS" w l ls 3, \
     "b_80_7.sim" notitle w l ls 2, \
     "b_80_7.test" notitle w l ls 1, \
     "SIMS/b80_7.dat" notitle w l ls 3, \
     "b_280_7.sim" notitle w l ls 2, \
     "b_280_7.test" notitle w l ls 1, \
     "SIMS/b280_7.dat" notitle w l ls 3, \
     "b_700_7.sim" notitle w l ls 2, \
     "b_700_7.test" notitle w l ls 1, \
     "SIMS/b700_7.dat" notitle w l ls 3, \
     "b_2400_7.sim" notitle w l ls 2, \
     "b_2400_7.test" notitle w l ls 1, \
     "SIMS/b2400_7.dat" notitle w l ls 3

