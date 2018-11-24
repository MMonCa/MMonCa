all: as_15_8.sim as_100_8.sim

as_15_8.sim: As_15_8.in
	$(IIS) $(PWD)/As_15_8.in
	cp DAs_15_8.in/As_15_8.01_01 as_15_8.sim
	rm -r DAs_15_8.in
	rm *txt

as_100_8.sim: As_100_8.in
	$(IIS) $(PWD)/As_100_8.in
	cp DAs_100_8.in/As_100_8.01_01 as_100_8.sim
	rm -r DAs_100_8.in
	rm *txt

plot: all As_Si_100_8.plot
	gnuplot As_Si_100_8.plot

clean:
	rm -rf as_15_8.sim as_100_8.sim D*.in Total_time

