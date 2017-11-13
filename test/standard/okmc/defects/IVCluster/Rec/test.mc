param set type=map<string,string>   key=MC/General/materials value="Silicon Si"

proc material { x y z } { return "Silicon" }

set size 20

param set type=bool             key=MC/Mesh/periodic.x value=true
param set type=map<string,bool> key=Silicon/Models/defined value=true index=IVCluster

init minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=IVCluster ID=VI   coord={ [expr $size/2+1] [expr $size/2]   [expr $size/2]   }
insert defect=IVCluster ID=V3I3 coord={ [expr $size/2]   [expr $size/2+1] [expr $size/2]   }
insert defect=IVCluster ID=V3I2 coord={ [expr $size/2]   [expr $size/2]   [expr $size/2+1] }

report all

test tag=b1 float=[extract count.particles defect=IVCluster ID=VI  ] value=2 error=0
test tag=b2 float=[extract count.particles defect=IVCluster ID=V3I3] value=6 error=0
test tag=b3 float=[extract count.particles defect=IVCluster ID=V3I2] value=5 error=0

anneal time=10 temp=500 events=1000000

test tag=a1 float=[extract count.particles particle=I] value=0 error=0
test tag=a2 float=[extract count.particles particle=V] value=1 error=0
test tag=a3 float=[extract count.particles] value=1 error=0
