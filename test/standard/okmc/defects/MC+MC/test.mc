param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas" 

proc tst { tag str1 str2 } {
     if { $str1 == $str2 } { 
     	lowmsg "$tag $str1 -- $str2 OK"
     } else {
       lowmsg "$tag $str1 != $str2 FAILED"
       error "$tag is $str1 instead of $str2 TEST FAILED"
     }
}

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

set size 15

param set type=proc      key=S_Iron/HeCluster/formation        value={ {
	set list "" 
        set eF 2.17
        append list "He2V2 < [expr 2*$eF -1.5] > "
	append list "He3V3 < [expr 3*$eF -6  ] > "
	append list "He4V4 < [expr 4*$eF -9  ] > "
	return $list
} }
param set type=proc     key=S_Iron/HeCluster/prefactor         value={ {
        set list "" 
        return $list
} }
param set type=proc     key=S_Iron/HeCluster/migration          value={ {
        set list "" 
        append list "He2V2 < 1e-3 1 > "
        append list "He3V3 < 0 5 > "
        append list "He4V4 < 0 5 > "
        return $list
} }
param set type=proc     key=S_Iron/HeCluster/transform.to          value={ {
        set list "" 
        append list "He2V2 < 0 5 > "
        append list "He3V3 < 0 5 > "
        append list "He4V4 < 0 5 > "
        return $list
} }
param set type=proc     key=S_Iron/HeCluster/transform.from          value={ {
        set list "" 
        append list "He2V2 < 0 5 > "
        append list "He3V3 < 0 5 > "
        append list "He4V4 < 0 5 > "
        return $list
} }



param set type=arrhenius             key=S_Iron/Vacancy/V(formation)            value={ 1e-5 10 }
param set type=arrhenius             key=S_Iron/Vacancy/V(migration)            value={ 1.0e-3 1.0 }
param set type=arrhenius             key=S_Iron/Iron/I(formation)        value={ 1e-5 10 }
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=He+HeV
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=HeV+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=HeCluster:He3V3+HeV 
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=V+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions        value=true      index=HeCluster+HeCluster
param set type=map<string,bool>      key=S_Iron/Models/defined             value=true      index=HeCluster

init minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

tst 1 [param get.reaction material=S_Iron index=HeCluster+HeCluster] true

insert defect=HeCluster ID=He2V2 coord={ [expr $size/2] [expr $size/2] [expr $size/2] }
insert defect=HeCluster ID=He2V2 coord={ [expr $size/2+1] [expr $size/2-1] [expr $size/2+1] }

test float=[extract count.defects defect=HeCluster ID=He2V2] value=2 error=0
test float=[extract count.defects defect=HeCluster ID=He4V4] value=0 error=0    

anneal time=10 temp=700 

test float=[extract count.defects defect=HeCluster ID=He2V2] value=0 error=0
test float=[extract count.defects defect=HeCluster ID=He4V4] value=1 error=0
