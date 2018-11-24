all: b_35_0.sim b_80_0.sim 

b_35_0.sim: B_35_0.in
	$(IIS) $(PWD)/B_35_0.in
	cp DB_35_0.in/B_35_0.01_01 b_35_0.sim
	rm -r DB_35_0.in
	rm *txt

b_80_0.sim: B_80_0.in
	$(IIS) $(PWD)/B_80_0.in
	cp DB_80_0.in/B_80_0.01_01 b_80_0.sim
	rm -r DB_80_0.in
	rm *txt

plot: all sse35.plot sse80.plot
	gnuplot sse35.plot
	gnuplot sse80.plot

clean:
	rm -rf b_35_0.sim b_80_0.sim D*.in  Total_time
