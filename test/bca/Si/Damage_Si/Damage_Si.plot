set term eps enh color

set output "Damage_Si.eps"

set xrange[0:900]
set yrange[1e16:1e22]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

set label 1 "15 keV, 8 10^{15} at/cm^{2}" at 50, 2e+21, 0 left norotate
set label 2 "80 keV, 8 10^{15} at/cm^{2}" at 400, 5e+20, 0 left norotate
set label 3 "10^{13} at/cm^{2}" at 40, 1.2e+18, 0 left norotate
set label 4 "80 keV, 10^{13} at/cm^{2}" at 300, 4e+17, 0 left norotate

plot "b_15_0_8e15.sim" t "Test" w l ls 2, \
     "SIMS/b15_0_8e15.dat" t "SIMS" w l ls 3, \
     "SIMS/b80_0_8e15.dat" notitle w l ls 3, \
     "b_80_0_8e15.sim" notitle w l ls 2, \
     "b_15_0.sim" notitle w l ls 2, \
     "b_80_0.sim" notitle w l ls 2, \
     "b_15_0_8e15.test" t "Ref" w l ls 1, \
     "b_80_0_8e15.test" notitle w l ls 1
