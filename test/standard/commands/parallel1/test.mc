param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"

param set type=int key=MC/General/domains    value=2
#param set type=int key=MC/General/subdomains value=2
#param set type=int key=MC/General/levels     value=3


#set time(600) 2.5e13
#set time(800) 6.5e8

set time(600) 1e12
set time(800) 6.5e7

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

set sizeX  8
set sizeYZ 80

param set type=arrhenius key=S_Iron/Iron/I(formation) value={ 1.5e25 3.7 }
param set type=arrhenius key=S_Iron/Iron/I(migration) value={ 1.0e-3 1.0 }
param set type=arrhenius key=S_Iron/Vacancy/V(formation)     value={ 1e-5 10 }
param set type=arrhenius key=S_Iron/Vacancy/V(migration)      value={ 5.0e-2 0.8 }
param set type=arrhenius key=S_Iron/Carbon/CI(migration)      value={ 1.0e-3   0.95 }
param set type=arrhenius key=S_Iron/Carbon/CI(formation)      value={ 175.51 2.83 }
param set type=array<string,string> key=S_Iron/Models/interactions value=true index=I+Gas


foreach { T } { 600 800 } {

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

report domains

set C_number 59
expr srand(555)
for { set n 0 } { $n < $C_number } { incr n } {
	insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=C
}

anneal time=$time($T) temp=$T

set kB  8.6174e-5
set If_P [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 0]
set If_E [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 1]
set CI   [expr $If_P*exp(-$If_E/($kB*($T+273.15)))]
test tag=$T.I array={ [extract profile.mobile name=I] } value=$CI init=1 end=$sizeX error=.10

set Im_P  [lindex [param get type=arrhenius key=S_Iron/Iron/I(migration)] 0]
set CC    [expr $C_number./$sizeX/$sizeYZ/$sizeYZ*1e21]

set Cif_P [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(formation)] 0]
set Cif_E [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(formation)] 1]
set Cf_P [lindex [param get type=arrhenius key=S_Iron/Carbon/C(formation)] 0]
set Cf_E [lindex [param get type=arrhenius key=S_Iron/Carbon/C(formation)] 1]
set CCi   [expr $CC*$Cif_P*exp(($Cf_E - $Cif_E)/($kB*($T+273.18)))]
set test  [test tag=$T.Ci array={ [extract profile.mobile name=CI] } value=$CCi init=1 end=$sizeX error=.25]

lowmsg "The averaged value is $test versus $CCi"

}

