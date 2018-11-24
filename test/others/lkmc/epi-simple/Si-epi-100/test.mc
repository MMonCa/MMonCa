param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true

param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=0.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=0.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride1 value=4e12   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride2 value=8e11   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride3 value=1e1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=2.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride1   value=2.4    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride2   value=1.9    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride3   value=0.1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.etch         value=0.0    index=Si

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

set vel(0)  51
set vel(9)  0.463
set vel(14) 3.332

set angle(0)  00
set angle(9)  55
set angle(14) 90

set sizeZ    [expr sqrt(2.)*.5431*26]
set sizeY    [expr sqrt(2.)*.5431*26]

param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true

set radians "$angle($i).0*2.0*3.1415926535897931/360.0"
set S [expr sin($radians)]
set C [expr cos($radians)]
set R [expr sqrt(2.0)]
set waferorient "i $C  j [expr $S/$R] k [expr $S/$R]"
set flatorient  "i -$S j [expr $C/$R] k [expr $C/$R]"

param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"

init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material

anneal time=1 temp=$T depth=51 epitaxy="Si 1.0"
set orig_time [extract time]
set orig_depth [lindex [extract ac.mean] 0]
lowmsg "Original depth is $orig_depth"

anneal time=1 temp=$T depth=31 epitaxy="Si 1.0"
set end_time [extract time]
set end_depth [lindex [extract ac.mean ] 0] 
lowmsg "$angle($i) - Final depth is $end_depth"
set velocity [expr ($orig_depth - $end_depth)/($end_time-$orig_time)*60]
set roughnes [extract ac.stdev]
lowmsg "Velocity  for angle $angle($i) is $velocity in nm/min."
lowmsg "Roughness for angle $angle($i) is $roughnes in nm."
test tag=$i float=$vel($i) value=$velocity error=0.03

save lammps=nodist-Si-100
