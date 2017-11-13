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

insert defect=<111> ID=I5 coord={ 5.5  5.5 5.5 }
insert defect=<111> ID=I5 coord={ 6.5  5.5 6.5 }
insert defect=<111> ID=I2 coord={ 3.5 5.5 5.5 }
insert defect=<111> ID=I2 coord={ 15.5 4.5 5.5 }
insert defect=<111> ID=I2 coord={ 5.5 5.5 7.5 }
insert defect=<111> ID=I3 coord={ 15.5 15.5 5.5 }
insert defect=ICluster ID=I4 coord={ 7 8 9 }

set compare "I2 3
I3 1
I5 2
"
test tag=val.0 one={[extract histogram defect=<111> material=S_Iron]} equal two={$compare}

set compare "I4 1
"
test tag=val.1 one={[extract histogram defect=ICluster material=S_Iron]} equal two={$compare}

set compare "I2 3
I3 1
I4 1
I5 2
"
test tag=val.2 one={[extract histogram defect=ICluster,<111> material=S_Iron]} equal two={$compare}
