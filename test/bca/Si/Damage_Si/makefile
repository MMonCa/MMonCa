all: b_15_0_8e15.sim b_80_0_8e15.sim

b_15_0_8e15.sim: B_15_0_8e15.in
	$(IIS) $(PWD)/B_15_0_8e15.in
	cp DB_15_0_8e15.in/B_15_0_8e15.01_01 b_15_0_8e15.sim
	rm -r DB_15_0_8e15.in
	rm *txt

b_80_0_8e15.sim: B_80_0_8e15.in
	$(IIS) $(PWD)/B_80_0_8e15.in
	cp DB_80_0_8e15.in/B_80_0_8e15.01_01 b_80_0_8e15.sim
	rm -r DB_80_0_8e15.in
	rm *txt

plot:	all Damage_Si.plot
	gnuplot Damage_Si.plot

clean:
	rm -rf b_15_0_8e15.sim b_80_0_8e15.sim D*.in Total_time
