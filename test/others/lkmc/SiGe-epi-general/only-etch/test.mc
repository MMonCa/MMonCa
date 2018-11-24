set T 600

proc material { x y z } { 
	if { ($x-12)*($x-12) + ($y -12)*($y-12) + ($z-12)*($z-12) < 2*2 } { 
		return "Silicon" }
	if { $x >14 } { return "Silicon" }
	return "Gas" 
}

#proc material { x y z } {
#	if { $x <18 && $y > 10 && $y < 14 } { return "Gas" }
#        if { $x >14 } { return "Silicon" }
#        return "Gas"
#} 

set sizeX [expr 24]
set sizeY [expr sqrt(2.)*0.5431*40]
set sizeZ [expr sqrt(2.)*0.5431*40]

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=bool                 key=Silicon/Epitaxy/model.simplified  value=false
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=1.0e9  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=1.0e6  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride1 value=4e12   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride2 value=8e11   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride3 value=1e1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=1.85    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride1   value=2.4    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride2   value=1.9    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride3   value=0.1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.etch         value=0.0    index=Si
param set type=proc                 key=Silicon/Epitaxy/prefactor.precursor value={ {
return 0 
} }

param set type=float                key=MC/General/snapshot.events        value=100000
param set type=int                  key=MC/General/snapshot.time.decade   value=1

proc snapshot { } { save lammps=nodist-Si append }

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material
save lammps=nodist-Si
anneal time=356000 temp=$T
save lammps=nodist-Si append

set end_depth [lindex [extract ac.mean] 0]
test tag=depth float=$end_depth value=17.602 error=0.02

set roughness [lindex [extract ac.stdev] 0]
test tag=rough float=$roughness value=0.09221 error=0.15
