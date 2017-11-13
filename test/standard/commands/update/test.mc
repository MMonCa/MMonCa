param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"

set time 1

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

set sizeX  8
set sizeYZ 25


param set type=arrhenius key=S_Iron/Vacancy/V(formation) value={ 3.0e25 3.9 }
param set type=arrhenius key=S_Iron/Iron/I(formation)    value={ 1.5e25 3.7 }
param set type=arrhenius key=S_Iron/Iron/I(migration)    value={ 1.0e-3 0.7 }
param set type=arrhenius key=S_Iron/Vacancy/V(migration) value={ 5.0e-2 0.8 }
param set type=map<string,float> key=Gas_S_Iron/Iron/recombination.length.right value={ I 1 }
param set type=map<string,float> key=Gas_S_Iron/Vacancy/recombination.length.right value={ V 1 }

param set type=int   key=MC/General/snapshot.time.decade value=10
param set type=float key=MC/General/snapshot.time.delta  value=10
param set type=float key=MC/General/snapshot.time.min    value=1e-5
param set type=float key=MC/General/snapshot.events value=1e6

set FILE [open "nodist-times.dat" "w"]
close $FILE

proc snapshot { } {
	set t [extract time]
	set FILE [open "nodist-times.dat" "a"]
	puts $FILE $t
	close $FILE
}

set T 1100

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

anneal time=1e25 temp=$T events=2e7

set kB  8.6174e-5

set VP [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 0]
set VE [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 1]
set VC [expr $VP*exp(-$VE/($kB*($T+273.15)))]
test tag=$T.V array={ [extract profile.mobile name=V] } value=$VC init=1 end=$sizeX error=.08

set IP [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 0]
set IE [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 1]
set IC [expr $IP*exp(-$IE/($kB*($T+273.15)))]
lowmsg "expr $IP*exp(-$IE/($kB*($T+273.15))) is $IC"
test tag=$T.I array={ [extract profile.mobile name=I] } value=$IC init=1 end=$sizeX error=.15

set rc [catch {exect diff nodist-times.dat Ref/times.dat} output]
if {$rc == 0 } {
 	test tag=two equal one=passed two=failed
} else {
	 test tag=two equal one=passed two=passed
}
