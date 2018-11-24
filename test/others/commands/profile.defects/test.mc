set size 50

proc material { x y z } { 
	return "Iron" 
}

param set type=map<string,string>   key=MC/General/materials value="Gas Gas Iron Fe"

#parameters
param set type=bool   key=MC/Mesh/periodic.x   value=true
param set type=float  key=MC/Mesh/spacing.x    value=1.3
param set type=float  key=MC/Mesh/delta.x value=1.
param set type=float  key=MC/Mesh/delta.y value=1.
param set type=float  key=MC/Mesh/delta.z value=1.

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I2 coord={ 5.5  5.5 5.5 }
insert defect=<100> ID=I2 coord={ 15.5 5.5 5.5 }

lowmsg [extract defects={ }]

lowmsg [extract profile name=Fe]

test tag=zero  array={ [extract profile name=Fe] } value=0 init=16 end=$size error=0
test tag=value array={ [extract profile name=Fe] } value=[expr 2./($size*$size*1e-21)] init=5 end=6 error=0      

test tag=val.1 array={ [extract profile name=Fe defect=<111>] } value=[expr 2./($size*$size*1e-21)] init=5 end=6 error=0
test tag=val.0 array={ [extract profile name=Fe defect=<100>] } value=0                             init=5 end=6 error=0
