set T 600


set value_depth 6.09
set value_rough 0.06

proc material { x y z } {
#	if { $x <7 && $y > 3 && $y < 5 } { return "Gas" }
        if { $x >6 } { return "Silicon" }
        return "Gas"
} 


set sizeX 50
set sizeY [expr sqrt(2.)*0.5431*8]
set sizeZ [expr sqrt(2.)*0.5431*8]

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=bool                 key=Silicon/Epitaxy/model.simplified  value=false
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch value=5.6e9  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig value=5  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=1.5    index=Si

param set type=float                key=MC/General/snapshot.events        value=500000
param set type=int                  key=MC/General/snapshot.time.decade   value=1


proc snapshot { } { save lammps=nodist-Si append no.print }
init minx=-20 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material

save lammps=nodist-Si
anneal time=5000 temp=$T epitaxy="Si 1.0"

save lammps=nodist-Si append

set rough [lindex [extract ac.stdev] 0]
set depth [lindex [extract no.print ac.mean] 0]
test tag=depth float=$depth value=$value_depth error=0.025
test tag=rough float=$rough value=$value_rough error=0.10
