param set type=map<string,string>   key=MC/General/materials value="GalliumArsenide GaAs Gas Gas"
param set type=map<string,string>   key=MC/Particles/elements value="As Arsenic,33,74.9216 Ga Gallium,11,22"
param set type=bool                 key=GalliumArsenide/Models/epitaxy value=true

param set type=int   key=MC/General/snapshot.time.decade value=1
param set type=float key=MC/General/snapshot.time.min    value=1
param set type=float key=MC/General/snapshot.events      value=250000

set T 600
set i 0 

set sizeZ    [expr sqrt(2.)*.5431*50]
set sizeY    [expr sqrt(2.)*.5431*50]


proc material { x y z } {
	global sizeY
	if { ($x-$sizeY/2)*($x-$sizeY/2) + ($y-$sizeY/2)*($y-$sizeY/2) + ($z-$sizeY/2)*($z-$sizeY/2) <10 } { return "GalliumArsenide" }
	return "Gas"
}


proc snapshot { } { save lammps=nodist-GaAs append }

init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material


save lammps=nodist-GaAs


anneal time=1e5 temp=$T events=10000000 epitaxy="Ga 1.0 As 1.0"

save lammps=nodist-GaAs append

set roughness [lindex [extract ac.stdev] 0]
test tag=rough value=$roughness float=5.16 error=.1
