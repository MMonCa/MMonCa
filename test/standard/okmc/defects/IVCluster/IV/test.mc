set size 20
param set type=map<string,string>   key=MC/General/materials value="Silicon Si"
param set type=float		key=Silicon/Models/amorphization.threshold value=-1
proc material { x y z } { return "Silicon" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=I coord={ 15 15 15 }
insert particle=V coord={ 15 15 15 }

report all

test tag=part.b float=[extract count.particles]                           value=2  error=0
test tag=defe.b float=[extract count.defects defect=IVCluster ID=VI] value=1  error=0

anneal time=1 temp=700

test tag=part.a float=[extract count.particles]                           value=0  error=0  
test tag=defe.a float=[extract count.defects defect=IVCluster ID=VI] value=0  error=0
