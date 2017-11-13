param set type=map<string,string>   key=MC/General/materials value="Gas Gas Silicon Si AmorphousSilicon aSi"

set T 550
set name extract

set i 0 

set amorph_x 50

proc material { x y z } {
     global amorph_x

        set res "Unknown"
	if { $x < 0 } { 
		set res "Gas" 
	} elseif { $x < $amorph_x } { 
		set res "AmorphousSilicon" 
	} else {
		set res "Silicon"
	}
	return $res
}

set sizeZ    [expr sqrt(2.)*.5431*30]
set sizeY    [expr sqrt(2.)*.5431*30]

param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true

set waferorient "i 1 j 0 k 0"
set flatorient  "i 0 j 1 k 1"

param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/flat.orientation  value="$flatorient"

init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material

set ac_mean   [extract ac.mean]
set ac_stdev  [extract ac.stdev]

test tag=ac.mean.x float=[lindex $ac_mean 0] value=$amorph_x error=0.05
test tag=ac.mean.y float=[lindex $ac_mean 1] value=[expr $sizeY / 2.] error=0.05
test tag=ac.mean.z float=[lindex $ac_mean 2] value=[expr $sizeZ / 2.] error=0.05

