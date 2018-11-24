set size 14.5 
set time 0.01
set temp 600

proc material { x y z } { return "SiliconCarbide" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true                            
param set type=map<string,string> key=MC/General/materials value={ SiliconCarbide SiC Gas Gas }
param set type=array<string,string> key=SiliconCarbide/Models/interactions value={ CI+SiC        IVCluster,1. } 
param set type=proc key=SiliconCarbide/IVCluster/prefactor value={ { 
set list ""
lappend list SiC^CI,CI     1.23e-3
return $list
} }
param set type=arrhenius key=SiliconCarbide/Silicon/SiI(migration) value={ 0 2 }
#param set type=arrhenius key=SiliconCarbide/Silicon/SiI(update) value={ 0 2 }

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material 

#inserting point defect
insert particle=SiI coord={7 7 7}

set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]
set Cl  [extract count.particles ID=SiC^CI   defect=ICluster]

test tag=SiC1 value=0 float=$SiC error=0
test tag=CSi1 value=0 float=$CSi error=0
test tag=CI1  value=0 float=$CI  error=0
test tag=VC1  value=0 float=$VC  error=0
test tag=SiI1 value=1 float=$SiI error=0
test tag=VSi1 value=0 float=$VSi error=0
test tag=Cl1  value=0 float=$Cl  error=0

anneal time=$time temp=$temp events=1


report all

set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI  defect=MobileParticle]
set VC  [extract count.particles particle=VC  defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]
set Cl  [extract count.particles ID=SiC^CI    defect=IVCluster]

test tag=SiC2 value=0 float=$SiC error=0
test tag=CSi2 value=0 float=$CSi error=0
test tag=CI2  value=0 float=$CI  error=0
test tag=VC2  value=0 float=$VC  error=0
test tag=SiI2 value=0 float=$SiI error=0
test tag=VSi2 value=0 float=$VSi error=0
test tag=Cl2  value=2 float=$Cl  error=0

anneal time=$time temp=$temp events=2

set Cl  [extract count.particles ID=SiC^CI    defect=IVCluster]

test tag=Cl3  value=0 float=$Cl  error=0
