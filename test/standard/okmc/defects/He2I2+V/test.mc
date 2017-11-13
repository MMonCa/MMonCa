param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

proc material { x y z } { return "S_Iron" }

set size 7

param set type=bool                  key=MC/Mesh/periodic.x value=true
param set type=arrhenius             key=S_Iron/Vacancy/V(formation)    value={ 1e-5 10 }
param set type=arrhenius             key=S_Iron/Vacancy/V(migration)    value={ 1.0e-3 1. }
param set type=arrhenius             key=S_Iron/Iron/I(migration)       value={ 1.0e-3 13.0 }
param set type=arrhenius             key=S_Iron/Iron/I(formation)       value={ 1e-5 10 }

param set type=array<string,string>  key=S_Iron/Models/interactions value=false index=V+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions value=true  index=HeCluster:HeI2+V
param set type=array<string,string>  key=S_Iron/Models/interactions value=true  index=HeCluster:He3I3+V

param set type=proc      key=S_Iron/HeCluster/formation           value={ {
	set list "" 
        set eF 10
        append list "HeI2  < [expr 2*$eF -.1] > "
	append list "He3I3 < [expr 3*$eF -.5] > "
	append list "He3I2 < [expr 3*$eF -.5] > "
	return $list
} }

param set type=proc     key=S_Iron/HeCluster/prefactor      value={ { return "HeI2,HeI 1e-3"   } }
param set type=proc     key=S_Iron/HeCluster/migration      value={ { return "HeI2     < 1e-3 1 > " } }
param set type=proc     key=S_Iron/HeCluster/transform.to   value={ { return "" } }
param set type=proc     key=S_Iron/HeCluster/transform.from value={ { return "" } }

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=HeCluster ID=HeI2 coord={ [expr $size/2]   [expr $size/2]   [expr $size/2] }
insert particle=V              coord={ [expr $size/2+1] [expr $size/2-1] [expr $size/2+1] }

test tag=HeI2.a float=[extract count.defects defect=Cluster ID=HeI2]            value=1 error=0
test tag=HeI2.a float=[extract count.defects defect=MobileParticle particle=V] value=1 error=0 
test tag=HeI2.a float=[extract count.defects defect=MobileParticle particle=HeI] value=0 error=0   

anneal time=10 temp=700 events=2e6

test tag=HeI2.b float=[extract count.defects defect=Cluster ID=HeI2]            value=0 error=0
test tag=HeI2.b float=[extract count.defects defect=MobileParticle particle=V] value=0 error=0 
test tag=HeI2.b float=[extract count.defects defect=MobileParticle particle=HeI] value=1 error=0   

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=HeCluster ID=He3I3 coord={ [expr $size/2] [expr $size/2] [expr $size/2] }
insert particle=V               coord={ [expr $size/2+1] [expr $size/2-1] [expr $size/2+1] }

test tag=He3I3.a float=[extract count.defects defect=Cluster ID=He3I3]           value=1 error=0
test tag=He3I3.a float=[extract count.defects defect=MobileParticle particle=V] value=1 error=0
test tag=He3I3.a float=[extract count.defects defect=Cluster ID=He3I2]           value=0 error=0

anneal time=10 temp=700 events=2e6

test tag=He3I3.b float=[extract count.defects defect=Cluster ID=He3I3]           value=0 error=0
test tag=He3I3.b float=[extract count.defects defect=MobileParticle particle=V] value=0 error=0
test tag=He3I3.b float=[extract count.defects defect=Cluster ID=He3I2]           value=1 error=0
