param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas AmorphousSilicon aSi"

proc material { x y z } { 
	if { $x < 13.5 } { return "Gas" }
	return "Silicon"
}

set T 550
set vel  3.332
set sizeZ  [expr sqrt(2.)*.5431*26]
set sizeY  90

param set type=bool key=MC/Mesh/periodic.y value=false
param set type=bool key=MC/Mesh/periodic.z value=true

set waferorient "i  0 j 1 k 1"
set flatorient  "i  1 j 0 k 0"
param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/flat.orientation  value="$flatorient"

param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+V

proc myName { x y z } {  if { $x < 52 } { return 3e22 } }

init minx=-3 miny=0 minz=0 maxx=60 maxy=$sizeY maxz=$sizeZ material=material

profile name=I		proc=myName

test array="[extract amorphous.fraction]" value=1 error=0 init=15 end=52  tag=AmorphIV 

anneal time=1 temp=$T depth=51
save lammps=nodist-SPER

set orig_time [extract time]
set orig_depth [lindex [extract ac.mean min.y=60 max.y=120 min.z=0 max.z=$sizeZ] 0]
lowmsg "Original depth is $orig_depth"

anneal time=1 temp=$T depth=31
set end_time [extract time]
set end_depth [lindex [extract ac.mean min.y=30 max.y=60 min.z=0 max.z=$sizeZ] 0] 
lowmsg "Final depth is $end_depth"
set velocity [expr ($orig_depth - $end_depth)/($end_time-$orig_time)*60]
set roughnes [extract ac.stdev]
lowmsg "Velocity  is $velocity in nm/min."
lowmsg "Roughness is $roughnes in nm."
test tag=vel float=$vel value=$velocity error=0.10

