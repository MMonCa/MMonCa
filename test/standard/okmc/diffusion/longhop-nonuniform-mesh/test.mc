param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas" 
set time 1000

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}


param set type=arrhenius key=S_Iron/Vacancy/V(formation) value={ 3.0e25 3.9 }
param set type=arrhenius key=S_Iron/Iron/I(formation)    value={ 1.5e25 3.7 }

param set type=arrhenius key=S_Iron/Iron/I(migration)    value={ 1.0e-3 1.0 }
param set type=arrhenius key=S_Iron/Vacancy/V(migration) value={ 5.0e-2 0.8 }
param set type=map<string,float> key=Gas_S_Iron/Iron/recombination.length.right value={ I 1 }
param set type=map<string,float> key=Gas_S_Iron/Vacancy/recombination.length.right value={ V 1 }
param set type=array<string,string> key=S_Iron/Models/interactions index=C+I value=false
param set type=array<string,string> key=S_Iron/Models/interactions index=C+V value=false

set T 1100

init linesx={-2 -1.5 -1 -0.5 0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0 6.4 6.8 7.2 7.6 8} linesy={0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25} linesz={0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25} material=material

#for { set x 0 } { $x < $sizeX } { incr x } {
#	for { set y 0 } { $y < $sizeYZ } { incr y } {
#		insert particle=C coord={$x $y 5}
#	}
#}


anneal time=$time temp=$T 

set kB  8.6174e-5

set VP [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 0]
set VE [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 1]
set VC [expr $VP*exp(-$VE/($kB*($T+273.15)))]
test tag=$T.V array={ [extract profile.mobile name=V] } value=$VC init=1 end=8 error=.01

set IP [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 0]
set IE [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 1]
set IC [expr $IP*exp(-$IE/($kB*($T+273.15)))]
test tag=$T.I array={ [extract profile.mobile name=I] } value=$IC init=1 end=8 error=.15

