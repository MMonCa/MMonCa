param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 
set time 1

proc material { x y z } {
	return "S_Iron"
}

set size 20

param set type=proc key=S_Iron/HeCluster/migration value={ { return "HeV2 < 1e-3 1 > " } }

param set type=array<string,string>  key=S_Iron/Models/interactions value=false      index=HeV+V
init minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
insert defect=HeCluster   ID=HeV2    coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

set one [extract defects={ }]
anneal time=1e5 temp=600 events=1000
set two [extract defects={ }]
test one={$one} not equal two={$two}

param set type=proc key=S_Iron/HeCluster/migration value={ { return "HeV2 < 0 1 > " } }

init minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
insert defect=HeCluster   ID=HeV2    coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

set one [extract defects={ }]
anneal time=1e5 temp=600 events=1000
set two [extract defects={ }]
test one={$one} equal two={$two}
