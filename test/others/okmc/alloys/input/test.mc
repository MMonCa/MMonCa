proc material { x y z } {
	if { $x < 0 } {
		return "Gas"
	}
	return "Iron"
}

set densityCm3 [param get type=float key=Iron/Models/alloy.density]
proc alloy { x y z } {
	global densityCm3
	if { 0 < $x && $x < 6 } {
		return 0
	} elseif { 6 < $x && $x < 12 } {
		return [expr 0.2 * $densityCm3]
	} elseif { 12 < $x && $x < 18 } {
		return [expr 0.4 * $densityCm3]
	} elseif { 18 < $x && $x < 24 } {
		return [expr 0.6 * $densityCm3]
	} elseif { 24 < $x && $x < 30 } {
		return [expr 0.8 * $densityCm3]
	} else {
		return $densityCm3
	}
}
param set type=map<string,string>   key=MC/General/materials value="Gas Gas Iron Fe"
param set type=array<string,string> key=Iron/Models/interactions index=V+V value=false
param set type=map<string,float>    key=Iron/Models/mixing.enthalpy value={ {xi 0 xo 0 x1 0 poly1 1 xi 1 xo 0} }

#param set type=arrhenius key=Iron/Iron/I(formation)    value={ 1 5 }
#param set type=arrhenius key=Iron/Iron/I(migration)    value={ 0% 9 0.9 5% 10 1.0 15% 11 1.1 75% 12 1.2 90% 13 1.3 100% 17 1.7 }
param set type=arrhenius key=Iron/Vacancy/V(migration) value={ 0.002 0.7 }
#param set type=arrhenius key=Iron/Helium/He(formation) value={ 0% 14 1.4 3% 15 1.5 98% 16 1.6 100% 18 1.8 }

set T 1100
set kB 8.6174e-5
proc rate { pref ener } {
	global T
	global kB
	return [expr $pref * exp ( - $ener / ($kB * ($T + 273.15)))]
}

# Formation Parameters test
set form0  { 5e26 3.0 }
set form02 { 4e26 3.0 }
set form04 { 3e26 3.0 }
set form06 { 2e26 3.0 }
set form08 { 1e26 3.0 }
set form1  { 9e25 2.9 }



set Vform "0% [lindex $form0 0] [lindex $form0 1]"
set Vform [concat $Vform "20% [lindex $form02 0] [lindex $form02 1]"]
set Vform [concat $Vform "40% [lindex $form04 0] [lindex $form04 1]"]
set Vform [concat $Vform "60% [lindex $form06 0] [lindex $form06 1]"]
set Vform [concat $Vform "80% [lindex $form08 0] [lindex $form08 1]"]
set Vform [concat $Vform "100% [lindex $form1 0] [lindex $form1 1]"]

param set type=arrhenius key=Iron/Vacancy/V(formation) value={ $Vform }

# Compute the equilibrium concentrations from rate theory
set C0  [rate [lindex $form0  0] [lindex $form0  1]]
set C02 [rate [lindex $form02 0] [lindex $form02 1]]
set C04 [rate [lindex $form04 0] [lindex $form04 1]]
set C06 [rate [lindex $form06 0] [lindex $form06 1]]
set C08 [rate [lindex $form08 0] [lindex $form08 1]]
set C1  [rate [lindex $form1  0] [lindex $form1  1]]

init minx=-1.5 miny=0 minz=0 maxx=36 maxy=12 maxz=12 material=material

profile name=Cr proc=alloy

anneal temp=$T time=1e25 events=5e7


test tag=Fe_0Cr array={ [extract profile.mobile name=V file=V.profile] } value=$C0 init=0 end=6 error=0.2
test tag=Fe_02Cr array={ [extract profile.mobile name=V] } value=$C02 init=6 end=12 error=0.2
test tag=Fe_04Cr array={ [extract profile.mobile name=V] } value=$C04 init=12 end=18 error=0.2
test tag=Fe_06Cr array={ [extract profile.mobile name=V] } value=$C06 init=18 end=24 error=0.2
test tag=Fe_08Cr array={ [extract profile.mobile name=V] } value=$C08 init=24 end=30 error=0.2
test tag=Fe_1Cr array={ [extract profile.mobile name=V] } value=$C1 init=30 end=36 error=0.2

