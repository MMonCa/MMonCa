param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 
set time 1

proc material { x y z } {
	return "S_Iron"
}

set size  8


param set type=array<string,string> key=S_Iron/Models/interactions value=true index=HeI+V
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=HeV+I

set Vf [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 1]
set Hef [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 1]
set HeVf [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)] 1]

lowmsg "Binding is [expr $Vf + $Hef - $HeVf]"

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=He coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

test tag=before.I error=0 value=1 float=[extract count.particles particle=He]
test tag=before.I error=0 value=0 float=[extract count.particles particle=HeI]
test tag=before.I error=0 value=0 float=[extract count.particles particle=V]

anneal time=1e20 events=1 temp=500

test tag=after.I error=0 value=0 float=[extract count.particles particle=He]
test tag=after.I error=0 value=1 float=[extract count.particles particle=HeI]
test tag=after.I error=0 value=1 float=[extract count.particles particle=V]


param set type=array<string,string> key=S_Iron/Models/interactions value=false index=HeI+V
param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeV+I

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=He coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

test tag=before.V error=0 value=1 float=[extract count.particles particle=He]
test tag=before.V error=0 value=0 float=[extract count.particles particle=HeV]
test tag=before.V error=0 value=0 float=[extract count.particles particle=I]

anneal time=1e20 events=1 temp=500

test tag=after.V error=0 value=0 float=[extract count.particles particle=He]
test tag=after.V error=0 value=1 float=[extract count.particles particle=HeV]
test tag=after.V error=0 value=1 float=[extract count.particles particle=I]

