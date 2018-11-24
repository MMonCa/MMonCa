set term eps enh color

set output "As_SiC.eps"

set xrange[0:500]
set yrange[1e16:5e19]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

set label 1 "40 keV" at 10, 3e+19, 0 left norotate
set label 2 "100 keV" at 80, 2e+19, 0 left norotate
set label 3 "300 keV" at 200, 1e+19, 0 left norotate
		
plot "as_300_16_30.sim" t "Test" w l ls 2, \
     "SIMS/As_SiC_300.sims" t "SIMS" w l ls 3, \
     "SIMS/As_SiC_40.sims" notitle w l ls 3, \
     "SIMS/As_SiC_100.sims" notitle w l ls 3, \
     "as_40_16_30.sim" notitle w l ls 2, \
     "as_sic_100kev.sim" notitle w l ls 2, \
     "as_300_16_30.test" t "Ref" w l ls 1, \
     "as_40_16_30.test" notitle w l ls 1, \
     "as_sic_100kev.test" notitle w l ls 1


