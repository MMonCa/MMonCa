param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas"

set T 550
set name Si-stress

proc material { x y z } {
        set res "Unknown"
	if { $x < 0 } { 
		set res "Gas" 
	} elseif { $x < 50 } { 
		set res "AmorphousSilicon" 
	} else {
		set res "Silicon"
	}
	return $res
}

set vel  [expr 0.7 * 10.43]

set sizeZ    [expr sqrt(2.)*.5431*40.]
set sizeY    [expr sqrt(2.)*.5431*40.]

param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true

set waferorient "i 1  j 0 k 0"
set flatorient  "i 0  j 1 k 1"

param set type=map<string,float> key=Silicon/Lattice/wafer.orientation          value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation           value="$flatorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/flat.orientation  value="$flatorient"

param set type=string key=Mechanics/General/model     value=Uniform
param set type=float  key=Mechanics/Uniform/stress.yy value=-1.e9

init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material

set orig_time  [extract time]
set orig_depth  [lindex [extract ac.mean] 0]
lowmsg "Original depth is $orig_depth"

anneal time=1 temp=$T depth=10

set end_time [extract time]
set end_depth  [lindex [extract ac.mean] 0] 

lowmsg "Final depth is $end_depth"

set velocity  [expr ($orig_depth - $end_depth)/($end_time-$orig_time)*60 ]
set roughness [extract ac.stdev]

lowmsg "Velocity  is $velocity nm/min."
lowmsg "Roughness is $roughness nm."

test tag=0 float=$vel value=$velocity error=0.10

