param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

set size 50

proc material { x y z } { 
	return "S_Iron" 
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x   value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I50 coord={ 49.25  25 49.25 }
extract defects={ } file=nodist-I50.dat
save lammps=nodist-I50

set rad1 [extract defect.radius defect=<111>]
test tag=radius.1 float=1.15233 value=$rad1 error=.0001

insert defect=<111> ID=I20 coord={ 15  15 15 }

set rad2 [extract defect.radius defect=<111> ID=I20]
test tag=radius.2 float=0.627822 value=$rad2 error=.0001

test tag=radius float=[expr ($rad1+$rad2)/2] value=[extract defect.radius defect=<111>] error=0.0001

test tag=min float=1.15233 value=[extract defect.radius defect=<111> min.radius=1] error=0.0001
