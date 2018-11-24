all: as_sic_100kev.sim as_40_16_30.sim as_300_16_30.sim

as_sic_100kev.sim: As_SiC_100keV.in
	$(IIS) $(PWD)/As_SiC_100keV.in
	cp DAs_SiC_100keV.in/As_SiC_100keV.01_01 as_sic_100kev.sim
	rm -r DAs_SiC_100keV.in
	rm *txt

as_40_16_30.sim: As_40_16_30.in
	$(IIS) $(PWD)/As_40_16_30.in
	cp DAs_40_16_30.in/As_40_16_30.01_01 as_40_16_30.sim
	rm -r DAs_40_16_30.in
	rm *txt 

as_300_16_30.sim: As_300_16_30.in
	$(IIS) $(PWD)/As_300_16_30.in
	cp DAs_300_16_30.in/As_300_16_30.01_01 as_300_16_30.sim
	rm -r DAs_300_16_30.in
	rm *txt 

plot :	all As_SiC.plot
	gnuplot As_SiC.plot

clean:
	rm -rf *.sim D*.in Total_time
