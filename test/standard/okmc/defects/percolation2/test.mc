set size 20
param set type=map<string,string>   key=MC/General/materials value="Silicon Si"
param set type=float		key=Silicon/Models/amorphization.threshold value=-1
proc material { x y z } { return "Silicon" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=I coord={ 15 15 15 }
insert particle=V coord={ 15 15 15 }

insert particle=I coord={ 15.6 15 15 }
insert particle=V coord={ 15.6 15 15 }

insert particle=V coord={ 15.3 15 15 }

report all

test tag=part.a float=[extract count.particles] value=5  error=0
test tag=defe.a float=[extract count.defects]   value=2  error=0

param set type=bool key=Silicon/IVCluster/percolation value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=I coord={ 15 15 15 }
insert particle=V coord={ 15 15 15 }

insert particle=I coord={ 15.6 15 15 }
insert particle=V coord={ 15.6 15 15 }

insert particle=V coord={ 15.3 15 15 }

report all

test tag=part.b float=[extract count.particles] value=5  error=0                     
test tag=defe.b float=[extract count.defects]   value=1  error=0
