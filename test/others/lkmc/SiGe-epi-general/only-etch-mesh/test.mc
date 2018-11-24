set T 600


param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=bool                 key=Silicon/Epitaxy/model.simplified  value=false
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=2.0e8  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=0.0  index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride1 value=4e12   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride2 value=8e11   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride3 value=1e1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=1.85    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride1   value=2.4    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride2   value=1.9    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride3   value=0.1    index=Si
param set type=proc                 key=Silicon/Epitaxy/prefactor.precursor  value={ {
return 0 
} }

param set type=float                key=MC/General/snapshot.events        value=5000
param set type=int                  key=MC/General/snapshot.time.decade   value=1

proc snapshot { } { save lammps=nodist-Si append }

init mesh=test.grd scale=30

save lammps=nodist-Si
anneal time=356000 temp=$T
save lammps=nodist-Si append


set end_depth [lindex [extract no.print ac.mean min.y=-20.0 max.y=-7.0 min.z=-20.0 max.z=14.0] 0]
set roughness [lindex [extract ac.stdev min.y=-20.0 max.y=-7.0 min.z=-20.0 max.z=14.0] 0]
lowmsg [format "%s %s %e %s %e" Substrate height $end_depth roughness $roughness]

test tag=depth float=$end_depth value=6.83 error=0.02
test tag=rough float=$roughness value=0.10 error=0.15
