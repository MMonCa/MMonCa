set T 700

param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=5.6e9  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=0  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=1.85    index=Si

proc material { x y z } { 
	if { ($x-12)*($x-12) + ($y -12)*($y-12) + ($z-12)*($z-12) < 10*10 } { 
		return "Silicon" }
	return "Gas" 
}

set sizeX [expr 40]
set sizeY [expr .3160*80]
set sizeZ [expr .3160*80]

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=int                  key=MC/General/random.seed value=451

#proc snapshot { } { save lammps=Si append no.print }
init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material
anneal time=1e15 temp=$T events=[expr 150000] epitaxy="Si 1."
save lammps=nodist-Si

set end_depth [lindex [extract ac.mean] 0]
set roughness [lindex [extract ac.stdev] 0]
lowmsg "depth: $end_depth roughness: $roughness"

test tag=depth float=$end_depth value=11.9688 error=0.03
test tag=rough float=$roughness value=3.8954 error=0.01
