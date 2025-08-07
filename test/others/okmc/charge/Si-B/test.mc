param set type=map<string,string>   key=MC/General/materials value="Silicon Si SiO2 SiO2"

set time 1

proc material { x y z } {
	if { $x < 0 } { return "SiO2" }
	return "Silicon"
}

set sizeX  12
set sizeYZ 40

param set type=bool      key=MC/Electrostatic/load.poisson   value=true
param set type=array<string,string> key=Silicon/Models/interactions value=true index=I+SiO2
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+BI
param set type=arrhenius key=Silicon/Vacancy/V(formation)    value="1 6"

set T 900
#set eF 0.33
set scale 0.720462
set kT [expr 8.6174e-5*($T+273.15)]

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

set B_number 59
expr srand(555)
for { set n 0 } { $n < $B_number } { incr n } {
	insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=B
}

anneal time=1 temp=$T events=1

set eF [expr -[extract electrostatic name=valence.bandedge dimension=0 min.x=0 max.x=$sizeX] ]
extract electrostatic file=nodist-cond.dat name=conduction.bandedge dimension=1 min.x=0 max.x=$sizeX
extract electrostatic file=nodist-elec.dat name=electron.density min.x=0 max.x=$sizeX
extract electrostatic file=nodist-hole.dat name=hole.density min.x=0 max.x=$sizeX
extract electrostatic file=nodist-band.dat name=bandgap min.x=0 max.x=$sizeX
extract electrostatic file=nodist-vale.dat name=valence.bandedge min.x=0 max.x=$sizeX

lowmsg "eF is $eF (versus 0.33)"

set kB  8.6174e-5
set If_P [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set If_E [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 2]
set CI   [expr $If_P*exp(-$If_E/($kB*($T+273.15)))]

set CB    [expr $B_number./$sizeX/$sizeYZ/$sizeYZ*1e21]

set Bif_P [lindex [param get type=arrhenius key=Silicon/Boron/BI(formation)] 0]
set Bif_E [lindex [param get type=arrhenius key=Silicon/Boron/BI(formation)] 1]
set Bf_P  [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 1]
set Bf_E  [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 2]

set eBiM  [param get type=float key=Silicon/Boron/BI(e(-1,0))]
set CBiM  [expr $CB*$Bif_P*exp(($Bf_E - $Bif_E)/($kB*($T+273.15)))]


anneal time=1e25 temp=$T events=2.4e8

test tag=$T.I array={ [extract profile.mobile name=I state=0] } value=$CI init=1 end=$sizeX error=.10
set test  [test tag=$T.BiM array={ [extract profile.mobile name=BI state=-] } value=$CBiM init=1 end=$sizeX error=.17]
lowmsg "The averaged value is $test versus $CBiM"

set CBi0 [expr $CBiM*exp(($scale*$eBiM - $eF)/($kB*($T+273.15)))]
test tag=$T.Bi0 array={ [extract profile.mobile name=BI state=0] } value=$CBi0 init=1 end=$sizeX error=.25
