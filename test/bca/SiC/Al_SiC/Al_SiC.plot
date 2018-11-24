set term eps enh color

set output "Al_SiC.eps"

set xrange[0:1400]
set yrange[1e16:4e19]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 lt 0 lc rgb "black"

set label 1 "30 keV" at 20, 2e+19, 0 left norotate
set label 2 "90 keV" at 50, 3e+19, 0 left norotate
set label 3 "195 keV" at 350, 1e+19, 0 left norotate
set label 4 "500 keV" at 500, 2e+18, 0 left norotate
set label 5 "1 MeV" at 900, 2e+18, 0 left norotate
		
plot "al_195_16_30.sim" t "Test" w l ls 2, \
     "al_30_16_30.sim" notitle w l ls 2, \
     "al_sic_90kev.sim" notitle w l ls 2, \
     "al_sic_500kev.sim" notitle w l ls 2, \
     "al_sic_1mev.sim" notitle w l ls 2, \
     "SIMS/Al_SiC_30.sims" t "SIMS" w l ls 3, \
     "SIMS/Al_SiC_195.sims" notitle w l ls 3, \
     "SIMS/Al_SiC_90keV.sims" notitle w l ls 3, \
     "SIMS/Al_SiC_1MeV.sims" notitle w l ls 3, \
     "SIMS/Al_SiC_500keV.sims" notitle w l ls 3, \
     "al_195_16_30.test" t "Ref" w l ls 1, \
     "al_30_16_30.test" notitle w l ls 1, \
     "al_sic_90kev.test" notitle w l ls 1, \
     "al_sic_500kev.test" notitle w l ls 1, \
     "al_sic_1mev.test" notitle w l ls 1

