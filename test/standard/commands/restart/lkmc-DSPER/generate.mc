param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas"

set T 550
set name SPER
#0 9 or 14
set i 0 
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
set vel(9)  0.463
set vel(14) 3.332

set angle(0)  00
set angle(9)  55
set angle(14) 90

set sizeZ    [expr sqrt(2.)*.5431*26]
set sizeY    [expr 60]

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

anneal time=1 temp=$T depth=50

restart save=dump

set time [extract time]
set depth [lindex [extract ac.mean] 0]
set roughness [lindex [extract ac.stdev] 0]

lowmsg "1.- $time $depth $roughness"

anneal time=5e3 temp=$T events=400000

set time [extract time]
set depth [lindex [extract ac.mean] 0]
set roughness [lindex [extract ac.stdev] 0]

lowmsg "2.- $time $depth $roughness"
