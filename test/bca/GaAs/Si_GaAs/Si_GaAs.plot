set term eps enh color

set output "Si_GaAs.eps"

set xrange[0:1000]
set yrange[1e16:3e18]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"

set label 1 "Si (7,30), 150 keV" at 50, 1e+17, 0 left norotate
set label 2 "Si (0,0), 150 keV" at 600, 2e+17, 0 left norotate

set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

plot "si150_gaasreo.sim" t "Test" w l ls 2, \
     "si150_gaas100.sim" notitle w l ls 2, \
     "SIMS/Si150_GaAsREO.dat" t "SIMS" w l ls 3, \
     "SIMS/Si150_GaAs100.dat"notitle w l ls 3, \
     "si150_gaasreo.test" t "Ref" w l ls 1, \
     "si150_gaas100.test" notitle w l ls 1

