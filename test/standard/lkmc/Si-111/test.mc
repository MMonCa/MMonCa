param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas"

set T 550
set name APL11
#0 9 or 14
set i 9
proc material { x y z } {
        set res "Unknown"
	if { $x < 0 } { 
		set res "Gas" 
	} elseif { $x < 52 } { 
		set res "AmorphousSilicon" 
	} else {
		set res "Silicon"
	}
	return $res
}

set vel(0)  10.43
set vel(9)  0.31
set vel(14) 3.332

set angle(0)  00
set angle(9)  55
set angle(14) 90

set sizeZ    [expr sqrt(2.)*.5431*26]
set sizeY    [expr 180]

param set type=bool key=MC/Mesh/periodic.y value=false
param set type=bool key=MC/Mesh/periodic.z value=true

set radians "$angle($i).0*2.0*3.1415926535897931/360.0"
set S [expr sin($radians)]
set C [expr cos($radians)]
set R [expr sqrt(2.0)]
set waferorient "i $C  j [expr $S/$R] k [expr $S/$R]"
set flatorient  "i -$S j [expr $C/$R] k [expr $C/$R]"

param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/flat.orientation  value="$flatorient"

init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material

anneal time=1 temp=$T depth=51
set orig_time [extract time]
set orig_depth  [lindex [extract ac.mean min.y=60 max.y=120 min.z=0 max.z=$sizeZ] 0]
lowmsg "Original depth is $orig_depth"

anneal time=1 temp=$T depth=31

save lammps=nodist-111

set end_time [extract time]
set end_depth [lindex [extract ac.mean min.y=60 max.y=120 min.z=0 max.z=$sizeZ] 0]
lowmsg "$angle($i) - Final depth is $end_depth"
set velocity [expr ($orig_depth - $end_depth)/($end_time-$orig_time)*60 ]
set roughnes [extract ac.stdev]
lowmsg "Velocity  for angle $angle($i) is $velocity in nm/min."
lowmsg "Roughness for angle $angle($i) is $roughnes in nm."
test tag=$i float=$vel($i) value=$velocity error=0.15

