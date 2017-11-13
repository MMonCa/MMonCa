param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

set size 20
set time 3000

proc material { x y z } { return "S_Iron" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=ICluster ID=I2 coord={ 5 5 5 }
insert defect=VCluster ID=V2 coord={ 15 15 15 }
insert defect=ICluster ID=I2 coord={ 15 5 5 }
insert defect=VCluster ID=V2 coord={ 5 15 15 }

report defects

test float=[extract count.particles]                   value=8 error=0
test float=[extract count.particles defect=Cluster]    value=8 error=0
test float=[extract count.particles defect=ICluster]   value=4 error=0
test float=[extract count.particles defect=VCluster]   value=4 error=0
test float=[extract count.positions position=I]        value=4 error=0
test float=[extract count.particles particle=V]                 value=4 error=0
test float=[extract count.positions position=I defect=ICluster] value=4 error=0
test float=[extract count.positions position=I defect=VCluster] value=0 error=0    

anneal time=$time temp=-81

test tag=1 float=[extract count.particles]               value=0 error=0

insert defect=ICluster ID=I2 coord={ 5 5 5 }
insert defect=VCluster ID=V3 coord={ 15 15 15 }
insert defect=ICluster ID=I2 coord={ 15 5 5 }
insert defect=VCluster ID=V2 coord={ 5 15 15 }

test tag=2.1 float=[extract count.particles]               value=9 error=0
test tag=2.2 float=[extract count.positions position=I]    value=4 error=0
test tag=2.3 float=[extract count.particles particle=V]    value=5 error=0

anneal time=$time temp=-81

report defects
test tag=3 float=[extract count.particles]               value=1 error=0
