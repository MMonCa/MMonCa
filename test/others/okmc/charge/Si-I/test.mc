param set type=map<string,string>   key=MC/General/materials value="Silicon Si SiO2 SiO2"

set time 1
set T 900
set eF 0.426088
set scale 0.720462
set kT [expr 8.6174e-5*($T+273.15)]

set eIM 0.7
set eIP 0.5
set eVM 0.8
set eVP 0.4

proc material { x y z } {
	if { $x < 0 } { return "SiO2" }
	return "Silicon"
}

set sizeX  10
set sizeYZ 25

param set type=bool      key=MC/Electrostatic/load.poisson value=true
param set type=float     key=Silicon/Silicon/I(e(-1,0))  value=$eIM
param set type=float     key=Silicon/Silicon/I(e(0,1))   value=$eIP
param set type=float     key=Silicon/Vacancy/V(e(-1,0))  value=$eVM
param set type=float     key=Silicon/Vacancy/V(e(0,1))   value=$eVP

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

set kB  8.6174e-5
set IP  [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set IE  [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 2]
set IC  [expr $IP*exp(-$IE/($kB*($T+273.15)))]
lowmsg "expr $IP*exp(-$IE/($kB*($T+273.15))) is $IC"
set ICM [expr $IC*exp(($eF-$eIM*$scale)/($kB*($T+273.15)))]
set ICP [expr $IC*exp(($eIP*$scale-$eF)/($kB*($T+273.15)))]

set VP [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 1]
set VE [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 2]
set VC [expr $VP*exp(-$VE/($kB*($T+273.15)))]
set VCM [expr $VC*exp(($eF-$eVM*$scale)/($kB*($T+273.15)))]
set VCP [expr $VC*exp(($eVP*$scale-$eF)/($kB*($T+273.15)))]

anneal time=1e25 temp=$T events=1e8

test tag=$T.I0 array={ [extract profile.mobile name=I state=0] } value=$IC  init=1 end=$sizeX error=.08
test tag=$T.I+ array={ [extract profile.mobile name=I state=+] } value=$ICP init=1 end=$sizeX error=.08
test tag=$T.I- array={ [extract profile.mobile name=I state=-] } value=$ICM init=1 end=$sizeX error=.08

test tag=$T.V0 array={ [extract profile.mobile name=V state=0] } value=$VC  init=1 end=$sizeX error=.08
test tag=$T.V+ array={ [extract profile.mobile name=V state=+] } value=$VCP init=1 end=$sizeX error=.08
test tag=$T.V- array={ [extract profile.mobile name=V state=-] } value=$VCM init=1 end=$sizeX error=.08
