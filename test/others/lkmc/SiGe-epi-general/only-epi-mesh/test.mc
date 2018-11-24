set T 600


set value_depthA 0.0
set value_depthB 6.1
set value_rough 0.0


#param set type=bool                key=MC/Mesh/periodic.x value=false
#param set type=bool                key=MC/Mesh/periodic.y value=false
#param set type=bool                key=MC/Mesh/periodic.z value=false

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=bool                 key=Silicon/Epitaxy/model.simplified  value=true
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.etch       value=0.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.mig        value=0.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride1 value=4e12   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride2 value=8e11   index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/prefactor.dehydride3 value=1e1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.precursor    value=2.0    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride1   value=2.4    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride2   value=1.9    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.dehydride3   value=0.1    index=Si
param set type=map<string,float>    key=Silicon/Epitaxy/barrier.etch         value=0.0    index=Si

param set type=float                key=MC/General/snapshot.events        value=5000
param set type=int                  key=MC/General/snapshot.time.decade   value=1


proc snapshot { } { save lammps=nodist-Si append no.print }
init mesh=../only-etch-mesh/test.grd scale=30
save lammps=nodist-Si

# Top of pillar
set roughA [lindex [extract ac.stdev min.y=-5.0 max.y=-1.0 min.z=-5.0 max.z=-1.0] 0]
set depthA [lindex [extract no.print ac.mean min.y=-5.0 max.y=-1.0 min.z=-5.0 max.z=-1.0] 0]
# Sample of substrate
set roughB [lindex [extract ac.stdev min.y=-20.0 max.y=-7.0 min.z=-20.0 max.z=14.0] 0]
set depthB [lindex [extract no.print ac.mean min.y=-18.0 max.y=-8.0 min.z=-18.0 max.z=12.0] 0]


test tag=depthA float=$depthA value=$value_depthA error=0.025
test tag=depthB float=$depthB value=$value_depthB error=0.025

extract profile name=Si dimension=2 file=test-Si2D.dat

lowmsg [format "%s %s %e %s %e" Pillar height $depthA roughness $roughA]
lowmsg [format "%s %s %e %s %e" Substrate height $depthB roughness $roughB]

anneal time=10.0 temp=$T epitaxy="Si 1.0"

save lammps=nodist-Si append



