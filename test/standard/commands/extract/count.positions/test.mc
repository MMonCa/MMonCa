param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"
set size 50

proc material { x y z } { 
	return "S_Iron" 
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x   value=true
param set type=float  key=MC/Mesh/spacing.x    value=1.3
param set type=float  key=MC/Mesh/delta.x value=1.
param set type=float  key=MC/Mesh/delta.y value=1.
param set type=float  key=MC/Mesh/delta.z value=1.

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111>    ID=I2 coord={ 5.5  5.5 5.5 }
insert defect=VCluster ID=V3 coord={ 15.5 5.5 5.5 }
insert particle=He     coord={ 5.5 5.5 15.5 }
insert particle=HeV    coord={ 5.5 15.5 5.5}

test tag=I  float=[extract count.positions position=I]  value=2 error=0
test tag=Fe float=[extract count.positions position=Fe] value=4 error=0      
test tag=V  float=[extract count.positions position=V]  value=1 error=0
