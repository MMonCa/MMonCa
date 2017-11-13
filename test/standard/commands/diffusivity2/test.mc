param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

set time 1

proc material { x y z } {
	return "S_Iron"
}

set sizeX  50
set sizeYZ 50

param set type=arrhenius key=S_Iron/Iron/I(formation)    value={ 1.5e25 3.7 }
param set type=arrhenius key=S_Iron/Iron/I(migration)    value={ 1.0e-3 1.0 }
param set type=arrhenius key=S_Iron/Vacancy/V(formation) value={ 1e-5 10 }
param set type=arrhenius key=S_Iron/Vacancy/V(migration) value={ 5.0e-2 0.8 }
param set type=arrhenius key=S_Iron/Carbon/CI(migration)   value={ 1.0e-3   0.95 }
param set type=arrhenius key=S_Iron/Carbon/CI(formation)   value={ 175.51 2.83 }
param set type=bool   key=MC/Mesh/periodic.x value=true
param set type=float  key=MC/Mesh/spacing.x value=2
param set type=float  key=MC/Mesh/spacing.y value=2
param set type=float  key=MC/Mesh/spacing.z value=2


foreach { T } { 600 } {

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

set C_number 59
expr srand(444)
for { set n 0 } { $n < $C_number } { incr n } {
	insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=C
}

insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=I

set kB  8.6174e-5

anneal time=1e25 temp=$T events=1e8

set CI   [extract profile.mobile name=I dimension=0]
lowmsg "Average I concentration is $CI"

set CC   [expr $C_number./$sizeX/$sizeYZ/$sizeYZ*1e21]
set If_P [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 0]
set If_E [lindex [param get type=arrhenius key=S_Iron/Iron/I(formation)] 1]
set CI0  [expr $If_P*exp(-$If_E/($kB*($T+273.15)))]

set Cif_P [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(formation)] 0]
set Cif_E [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(formation)] 1]
set Cim_P [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(migration)] 0]
set Cim_E [lindex [param get type=arrhenius key=S_Iron/Carbon/CI(migration)] 1]
set Cf_P [lindex [param get type=arrhenius key=S_Iron/Carbon/C(formation)] 0]
set Cf_E [lindex [param get type=arrhenius key=S_Iron/Carbon/C(formation)] 1]

set CCi   [expr $CI/$CI0*$CC*$Cif_P*exp(($Cf_E - $Cif_E)/($kB*($T+273.15)))]
set thDiff [expr $Cim_P * exp(-$Cim_E/($kB*($T+273.15))) * $CCi/($CC)]

set test  [test tag=$T.Ci array={ [extract profile.mobile name=CI dimension=1] } value=$CCi init=1 end=$sizeX error=.25]
test tag=$T.diff float=[extract diffusivity macroscopic material=S_Iron name=C] value=$thDiff error=.15

}
