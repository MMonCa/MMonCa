param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"
param set type=bool key=MC/General/average.rates value=true

set size 50
set fluence 1e12

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false

proc snapshots { } {
 report events
}

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

cascade file=cascade periodic fluence=$fluence flux=1e12 temp=-150 voluminic
report all
test tag=standard float=[extract count.particles] value=[expr $fluence*$size*$size*1e-14] error=0

