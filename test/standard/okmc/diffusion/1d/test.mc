param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

set time 1

proc material { x y z } {
	return "S_Iron"
}

set sizeX  50
set sizeYZ 50

param set type=bool    key=MC/Mesh/periodic.x value=true
param set type=float   key=MC/Mesh/spacing.x value=2
param set type=float   key=MC/Mesh/spacing.y value=2
param set type=float   key=MC/Mesh/spacing.z value=2
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=I+I
param set type=coordinates key=S_Iron/Iron/I(axis) value="0 0 0" new

set I_number 3
set CI   [expr $I_number./$sizeX/$sizeYZ/$sizeYZ*1e21]

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material


expr srand(555)
for { set n 0 } { $n < $I_number } { incr n } {
	insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=I
}

anneal time=1e25 temp=600 events=5e7

set test  [test tag=3D array={ [extract profile.mobile name=I] } value=$CI init=1 end=$sizeX error=.075]

#---------------------------
param set type=coordinates key=S_Iron/Iron/I(axis) value="1 0 0"

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

expr srand(555)
for { set n 0 } { $n < $I_number } { incr n } {
	insert coord={ [expr rand()*$sizeX] [expr rand()*$sizeYZ] [expr rand()*$sizeYZ] } particle=I
}

anneal time=1e25 temp=600 events=5e7

set test  [test tag=1D array={ [extract profile.mobile name=I] } not value=$CI init=1 end=$sizeX error=2]

lowmsg [extract defects={ }]
