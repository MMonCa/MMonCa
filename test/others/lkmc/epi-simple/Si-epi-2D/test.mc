param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=int                  key=MC/General/random.seed value=56478
set value(0) 35.3397 
set value(5) 33.1694 
set value(10) 32.7834 
set value(15) 32.6078 
set value(20) 32.2719 
set value(25) 32.9806 
set value(30) 32.3103 
set value(35) 32.66 
set value(40) 32.3055 
set value(45) 32.5393 
set value(50) 32.3302 
set value(55) 32.5487 
set value(60) 32.5494 
set value(65) 32.7547 
set value(70) 32.5738 
set value(75) 32.3406 
set value(80) 32.3533 
set value(85) 33.1842 
set value(90) 33.4845 
set value(95) 35.0307


set T 600
set name SPER
#0 9 or 14
set i 0 
proc material { x y z } {
	if { $x < 52 } { 
		set res "Gas" 
	} else {
		set res "Silicon"
	}
	return $res
}

set vel(0)  5
set vel(9)  0.463
set vel(14) 3.332

set angle(0)  00
set angle(9)  55
set angle(14) 90

set sizeZ    [expr sqrt(2.)*.5431*26]
set sizeY    100

param set type=bool key=MC/Mesh/periodic.y value=false
param set type=bool key=MC/Mesh/periodic.z value=true

set radians "$angle($i).0*2.0*3.1415926535897931/360.0"
set S [expr sin($radians)]
set C [expr cos($radians)]
set R [expr sqrt(2.0)]
set waferorient "i $C  j [expr $S/$R] k [expr $S/$R]"
set flatorient  "i -$S j [expr $C/$R] k [expr $C/$R]"

param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=0.0  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=0.0  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.hydride    value=0.0   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride1 value=4e12   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride2 value=8e11   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride3 value=1e1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=2.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride1   value=2.4    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride2   value=1.9    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride3   value=0.1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.etch         value=0.0    index=Si


param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"


init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material

anneal time=1 temp=$T depth=31 epitaxy="Si 1."

save lammps=nodist-Si-epi-2D

set roughness [lindex [extract ac.stdev] 0]

foreach i "0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95" {
        set depth [lindex [extract no.print ac.mean min.y=$i max.y=[expr $i+5]] 0]
#	lowmsg "set value($i) $depth"
        test tag=depth.$i float=$depth value=$value($i)0 error=0.01
}

test tag=rough float=$roughness value=2.4 error=0.02
