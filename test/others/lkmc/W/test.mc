set T 25

proc material { x y z } { 
	if { ($x-12)*($x-12) + ($y -12)*($y-12) + ($z-12)*($z-12) < 2*2 } { 
		return "Tungsten" }
	return "Gas" 
}

set sizeX [expr 24]
set sizeY [expr .3160*80]
set sizeZ [expr .3160*80]

param set type=map<string,string>   key=MC/General/materials value="Tungsten W Gas Gas"
param set type=map<string,string>   key=MC/Particles/elements value={ W Tungsten,74,183.84 }
param set type=map<string,bool>     key=Tungsten/Models/particles value={ }
param set type=map<string,bool>     key=Tungsten/Models/defined value={ }
param set type=array<string,string> key=Tungsten/Models/interactions value={ }
param set type=array<string>        key=Tungsten/Models/interaction.result value={ } new
param set type=bool                 key=Tungsten/Models/epitaxy value=true

set waferorient "i 0 j 0 k 1"
set flatorient  "i 0 j 1 k 0"
param set type=map<string,float> key=Tungsten/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Tungsten/Lattice/flat.orientation    value="$flatorient"

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material

lowmsg [extract coordination coord={ 12 12 12 } radius=.6 type=1]
anneal time=1e15 temp=$T events=[expr 150000] epitaxy="W 1".0
save lammps=nodist-W-end

set roughness [lindex [extract ac.stdev] 0]
set depth [lindex [extract no.print ac.mean] 0]

test tag=depth float=$depth     value=11.9606 error=0.008
test tag=rough float=$roughness value=4.30556 error=0.005

