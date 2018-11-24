all: al_30_16_30.sim al_195_16_30.sim al_sic_40kev.sim al_sic_90kev.sim al_sic_500kev.sim al_sic_1mev.sim

al_30_16_30.sim: Al_30_16_30.in
	$(IIS) $(PWD)/Al_30_16_30.in
	cp DAl_30_16_30.in/Al_30_16_30.01_01 al_30_16_30.sim
	rm -r DAl_30_16_30.in
	rm *txt

al_195_16_30.sim: Al_195_16_30.in
	$(IIS) $(PWD)/Al_195_16_30.in
	cp DAl_195_16_30.in/Al_195_16_30.01_01 al_195_16_30.sim
	rm -r DAl_195_16_30.in
	rm *txt

al_sic_40kev.sim: Al_SiC_40keV.in
	$(IIS) $(PWD)/Al_SiC_40keV.in
	cp DAl_SiC_40keV.in/Al_SiC_40keV.01_01 al_sic_40kev.sim
	rm -r DAl_SiC_40keV.in
	rm *txt

al_sic_90kev.sim: Al_SiC_90keV.in
	$(IIS) $(PWD)/Al_SiC_90keV.in
	cp DAl_SiC_90keV.in/Al_SiC_90keV.01_01 al_sic_90kev.sim
	rm -r DAl_SiC_90keV.in
	rm *txt

al_sic_500kev.sim: Al_SiC_500keV.in
	$(IIS) $(PWD)/Al_SiC_500keV.in
	cp DAl_SiC_500keV.in/Al_SiC_500keV.01_01 al_sic_500kev.sim
	rm -r DAl_SiC_500keV.in
	rm *txt

al_sic_1mev.sim: Al_SiC_1MeV.in
	$(IIS) $(PWD)/Al_SiC_1MeV.in
	cp DAl_SiC_1MeV.in/Al_SiC_1MeV.01_01 al_sic_1mev.sim
	rm -r DAl_SiC_1MeV.in
	rm *txt

plot:	all	Al_SiC.plot
	gnuplot Al_SiC.plot

clean:
	rm -rf *.sim D*.in Total_time
