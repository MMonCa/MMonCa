param set type=map<string,string>   key=MC/General/materials value="GalliumArsenide GaAs Gas Gas"
param set type=map<string,string>   key=MC/Particles/elements value="As Arsenic,33,74.9216 Ga Gallium,11,22"
param set type=array<string,string> key=GalliumArsenide/Models/interactions value={ }
param set type=map<string,bool>     key=GalliumArsenide/Models/particles value={As false}
param set type=bool                 key=GalliumArsenide/Models/epitaxy value=true

param set type=int   key=MC/General/snapshot.time.decade value=1
param set type=float key=MC/General/snapshot.time.min    value=100
param set type=float key=MC/General/snapshot.events      value=2000000

set T 600
set i 0 

proc material { x y z } {
	if { $x > 52 } { return "GalliumArsenide" }
	if { $x > 47 && $y > 10 && $y < 25 && $z > 10 && $z < 14 } {
		return "GalliumArsenide" } 
	return "Gas"
}

set sizeZ    [expr sqrt(2.)*.5431*50]
set sizeY    [expr sqrt(2.)*.5431*50]

param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true

set waferorient "i 1 j 0 k 0"
set flatorient  "i 0 j 1 k 1"

param set type=map<string,float> key=GalliumArsenide/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=GalliumArsenide/Lattice/flat.orientation  value="$flatorient"

proc snapshot { } { save lammps=nodist-GaAs-100 append }

init minx=-2 miny=0 minz=0 maxx=58 maxy=$sizeY maxz=$sizeZ material=material

save lammps=nodist-GaAs-100


anneal time=1 temp=$T depth=31 events=40000000 epitaxy="Ga 1.0 As 1.0"

save lammps=nodist-GaAs-100 append


set roughness  [extract ac.stdev]

lowmsg $roughness

test tag=rough0 value=[lindex $roughness 0] float=6.67 error=.1
test tag=rough1 value=[lindex $roughness 1] float=10.6 error=.1
test tag=rough2 value=[lindex $roughness 2] float=11.1 error=.1

