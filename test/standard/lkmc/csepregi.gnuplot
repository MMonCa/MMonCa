!grep elocity Si-1*/*.log > SPER.txt

set term postscript eps enh color 22
set output "csepregi.eps"

set size 0.68,0.60

set yrange [-0.5:12.5]
set xtics ("0^o" 0, "10^o" 10, "20^o" 20, "30^o" 30, "40^o" 40, "50^o" 50, "60^o" 60, "70^o" 70, "80^o" 80, "90^o" 90)
set xlabel "{/Symbol q} (degrees)"
set ylabel "Growth velocity v (nm/min)" offset 2,-1

set label "(100)" at  0,4.        font "Times-Roman,18"
set label "(311)" at 25,4. center font "Times-Roman,18"
set label "(111)" at 55,4. center font "Times-Roman,18"
set label "(011)" at 81,4.        font "Times-Roman,18"

set pointsize 1.5
set tmargin .4
set rmargin .4

plot \
"csepregi.dat"     w l lt 1 lw 1 t "Csepregi et al.", \
"SPER.txt" u 4:6       w p lt 3 lw 10 t "LKMC"

