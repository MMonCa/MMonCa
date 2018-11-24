set term eps enh color

set output "Se_GaAs.eps"

set xrange[0:1000]
set yrange[1e16:3e18]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"

set label 1 "Se (7,30), 300 keV" at 250, 2e+18, 0 left norotate
set label 2 "Se (0,0), 300 keV" at 400, 2e+17, 0 left norotate

set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

plot "se300_gaasreo.sim" t "Test" w l ls 2, \
     "se300_gaas100.sim" notitle w l ls 2, \
     "SIMS/Se300_GaAsREO.dat" t "SIMS" w l ls 3, \
     "SIMS/Se300_GaAs100.dat" notitle w l ls 3, \
     "se300_gaasreo.test" t "Ref" w l ls 1, \
     "se300_gaas100.test" notitle w l ls 1

