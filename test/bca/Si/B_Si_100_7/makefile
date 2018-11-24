all: b_15_7.sim b_80_7.sim b_280_7.sim b_700_7.sim b_2400_7.sim

b_15_7.sim: B_15_7.in
	$(IIS) $(PWD)/B_15_7.in
	cp DB_15_7.in/B_15_7.01_01 b_15_7.sim
	rm -r DB_15_7.in
	rm *txt 

b_80_7.sim: B_80_7.in
	$(IIS) $(PWD)/B_80_7.in
	cp DB_80_7.in/B_80_7.01_01 b_80_7.sim
	rm -r DB_80_7.in
	rm *txt 

b_280_7.sim: B_280_7.in
	$(IIS) $(PWD)/B_280_7.in
	cp DB_280_7.in/B_280_7.01_01 b_280_7.sim
	rm -r DB_280_7.in
	rm *txt 

b_700_7.sim: B_700_7.in
	$(IIS) $(PWD)/B_700_7.in
	cp DB_700_7.in/B_700_7.01_01 b_700_7.sim
	rm -r DB_700_7.in
	rm *txt 

b_2400_7.sim: B_2400_7.in
	$(IIS) $(PWD)/B_2400_7.in
	cp DB_2400_7.in/B_2400_7.01_01 b_2400_7.sim
	rm -r DB_2400_7.in
	rm *txt 

plot: all B_Si_100_7.plot
	gnuplot B_Si_100_7.plot

clean:
	rm -rf *.sim D*.in Total_time
