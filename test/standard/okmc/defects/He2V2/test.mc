param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

proc material { x y z } {
	return "S_Iron"
}

proc tst { tag str1 str2 } {
     if { $str1 == $str2 } {
        lowmsg "$tag $str1 -- $str2 OK"
     } else {
       lowmsg "$tag $str1 != $str2 FAILED"
       error "$tag is $str1 instead of $str2 TEST FAILED"
     }
}

set size 8

param set type=bool                  key=MC/Mesh/periodic.x        value=trye
param set type=arrhenius             key=S_Iron/Vacancy/V(formation) value={ 1e-5 5 }
param set type=arrhenius             key=S_Iron/Vacancy/V(migration) value={ 1.0e-3 1.0 }
param set type=arrhenius             key=S_Iron/Iron/I(formation)    value={ 5e-5 10 }
param set type=array<string,string>  key=S_Iron/Models/interactions  value=false   index=He+HeV
param set type=array<string,string>  key=S_Iron/Models/interactions  value=false   index=HeV+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions  value=false   index=V+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions  value=true    index=HeCluster:He2V+V
param set type=map<string,bool>      key=S_Iron/Models/defined       value=true    index=HeCluster
param set type=proc                  key=S_Iron/HeCluster/prefactor  value={ }

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

tst a [param get.reaction material=S_Iron index=HeCluster:He2V+V] true

insert defect=HeCluster ID=He2V    coord={ [expr $size/2] [expr $size/2] [expr $size/2] }
insert                  particle=V coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

anneal time=.1 temp=700 

test float=[extract count.particles defect=Cluster ID=He2V2] value=4 error=0
