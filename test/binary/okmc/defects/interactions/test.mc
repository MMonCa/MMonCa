set size 14.5 
set time 0.01
set temp 600

proc material { x y z } { return "SiliconCarbide" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true                            
param set type=map<string,string> key=MC/General/materials value={ SiliconCarbide SiC Gas Gas }
param set type=array<string,string> key=SiliconCarbide/Models/interactions value={ {
  SiI+CSi       true
  SiI+VSi       0   
  CI+VC         0   
  SiI+VC        true
  CI+VSi        true  
  CI+SiC        true
  SiC+VSi       true
} }

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material 

#inserting point defect
insert particle=VSi coord={7 7 7}
insert particle=SiC coord={7.1 7.1 7.1}
anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

test tag=SiC1 value=0 float=$SiC error=0
test tag=CSi1 value=0 float=$CSi error=0
test tag=CI1  value=0 float=$CI  error=0
test tag=VC1  value=1 float=$VC  error=0
test tag=SiI1 value=0 float=$SiI error=0
test tag=VSi1 value=0 float=$VSi error=0

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

#inserting point defect
insert particle=VSi coord={7 7 7}
insert particle=CI coord={7.1 7.1 7.1}

report all

anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI  defect=MobileParticle]
set VC  [extract count.particles particle=VC  defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

report defects reactions

test tag=SiC2 value=0 float=$SiC error=0
test tag=CSi2 value=1 float=$CSi error=0
test tag=CI2  value=0 float=$CI  error=0
test tag=VC2  value=0 float=$VC  error=0
test tag=SiI2 value=0 float=$SiI error=0
test tag=VSi2 value=0 float=$VSi error=0

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

#inserting point defect
insert particle=VSi coord={7 7 7}
insert particle=SiI coord={7.1 7.1 7.1}
anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

test tag=SiC3 value=0 float=$SiC error=0
test tag=CSi3 value=0 float=$CSi error=0
test tag=CI3  value=0 float=$CI  error=0
test tag=VC3  value=0 float=$VC  error=0
test tag=SiI3 value=0 float=$SiI error=0
test tag=VSi3 value=0 float=$VSi error=0

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

#inserting point defect
insert particle=VC coord={7 7 7}
insert particle=CI coord={7.1 7.1 7.1}
anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

test tag=SiC4 value=0 float=$SiC error=0
test tag=CSi4 value=0 float=$CSi error=0
test tag=CI4  value=0 float=$CI  error=0
test tag=VC4  value=0 float=$VC  error=0
test tag=SiI4 value=0 float=$SiI error=0
test tag=VSi4 value=0 float=$VSi error=0

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

#inserting point defect
insert particle=SiC coord={7 7 7}
insert particle=CI coord={7.1 7.1 7.1}
anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

test tag=SiC5 value=0 float=$SiC error=0
test tag=CSi5 value=0 float=$CSi error=0
test tag=CI5  value=0 float=$CI  error=0
test tag=VC5  value=0 float=$VC  error=0
test tag=SiI5 value=1 float=$SiI error=0
test tag=VSi5 value=0 float=$VSi error=0

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

#inserting point defect
insert particle=SiI coord={7 7 7}
insert particle=CSi coord={7.1 7.1 7.1}
anneal time=$time temp=$temp
set SiC [extract count.particles particle=SiC defect=MobileParticle]
set CSi [extract count.particles particle=CSi defect=MobileParticle]
set CI  [extract count.particles particle=CI defect=MobileParticle]
set VC  [extract count.particles particle=VC defect=MobileParticle]
set SiI [extract count.particles particle=SiI defect=MobileParticle]
set VSi [extract count.particles particle=VSi defect=MobileParticle]

test tag=SiC6 value=0 float=$SiC error=0
test tag=CSi6 value=0 float=$CSi error=0
test tag=CI6  value=1 float=$CI  error=0
test tag=VC6  value=0 float=$VC  error=0
test tag=SiI6 value=0 float=$SiI error=0
test tag=VSi6 value=0 float=$VSi error=0

