set term eps enh color

set output "sse80.eps"

set xrange[0:1400]
set yrange[1e15:1e18]
set logscale y
set format y "10^{%T}"
set xlabel "Depth [nm]"
set ylabel "Concentration [cm^{-3}]"
set style line 1 lt 1 lw 3 lc rgb "black"
set style line 2 lt 1 lc rgb "red"
set style line 3 pt 6 lc rgb "black"
set style line 4 lt 0 lc rgb "black" 
set style line 5 lt 1 lc rgb "black"

plot "Data/b80_0.dat" smooth csplines t "SIMS" w p ls 3, \
     "Data/B80_0.out" t "MARLOWE" w l ls 4, \
     "Data/B80_0.cai" t "UT-MARLOWE" w l ls 5, \
     "b_80_0.sim" t "Test" w l ls 2, \
     "b_80_0.test" t "Ref" w l ls 1 



