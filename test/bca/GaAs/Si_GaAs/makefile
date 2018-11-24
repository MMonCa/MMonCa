all: si150_gaasreo.sim si150_gaas100.sim

si150_gaasreo.sim: Si150_GaAsREO.in
	$(IIS) $(PWD)/Si150_GaAsREO.in
	cp DSi150_GaAsREO.in/Si150_GaAsREO.01_01 si150_gaasreo.sim
	rm -r DSi150_GaAsREO.in
	rm *txt 

si150_gaas100.sim: Si150_GaAs100.in
	$(IIS) $(PWD)/Si150_GaAs100.in
	cp DSi150_GaAs100.in/Si150_GaAs100.01_01 si150_gaas100.sim
	rm -r DSi150_GaAs100.in
	rm *txt 

plot:	all	Si_GaAs.plot
	gnuplot Si_GaAs.plot

clean:
	rm -rf *.sim D*.in Total_time
