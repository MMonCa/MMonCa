all: se300_gaasreo.sim se300_gaas110.sim se300_gaas100.sim

se300_gaasreo.sim: Se300_GaAsREO.in
	$(IIS) $(PWD)/Se300_GaAsREO.in
	cp DSe300_GaAsREO.in/Se300_GaAsREO.01_01 se300_gaasreo.sim
	rm -r DSe300_GaAsREO.in
	rm *txt 

se300_gaas110.sim: Se300_GaAs110.in
	$(IIS) $(PWD)/Se300_GaAs110.in
	cp DSe300_GaAs110.in/Se300_GaAs110.01_01 se300_gaas110.sim
	rm -r DSe300_GaAs110.in
	rm *txt 

se300_gaas100.sim: Se300_GaAs100.in
	$(IIS) $(PWD)/Se300_GaAs100.in
	cp DSe300_GaAs100.in/Se300_GaAs100.01_01 se300_gaas100.sim
	rm -r DSe300_GaAs100.in
	rm *txt 

plot:	all	Se_GaAs.plot
	gnuplot Se_GaAs.plot

clean:
	rm -rf *.sim D*.in Total_time
