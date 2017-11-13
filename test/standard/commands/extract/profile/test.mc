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

insert particle=FeI coord={ 5.5 5.5 5.5 }

test tag=zero  array={ [extract profile name=FeI] } value=0 init=7 end=$size error=0
test tag=value array={ [extract profile name=FeI] } value=[expr 1./($size*$size*1e-21)] init=5 end=6 error=0      

test tag=val.2 array.2D={ [extract profile name=FeI dimension=2] } value=[expr 1./($size*1e-21)] init={ 5 5 5 } end={ 6 6 6 } error=0
test tag=val.3 array.3D={ [extract profile name=FeI dimension=3] } value=[expr 1./(1e-21)] init={ 5 5 5 } end={ 6 6 6 } error=0

insert particle=FeI coord={ 5.7 5.3 5.2 }

test tag=val.10    array={ [extract profile name=FeI dimension=1] } value=[expr 2./($size*$size*1e-21)] init=5 end=6 error=0
test tag=val.20 array.2D={ [extract profile name=FeI dimension=2] } value=[expr 2./($size*1e-21)] init={ 5 5 5 } end={ 6 6 6 } error=0
test tag=val.30 array.3D={ [extract profile name=FeI dimension=3] } value=[expr 2./(1e-21)] init={ 5 5 5 } end={ 6 6 6 } error=0
