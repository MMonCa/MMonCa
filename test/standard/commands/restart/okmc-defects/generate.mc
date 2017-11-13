param set type=map<string,string>   key=MC/General/materials value="Silicon Si SiO2 SiO2"
param set type=float		    key=Silicon/Models/amorphization.threshold value=-1
set time 1

proc material { x y z } { 
	if { $x < 12 } { return "Silicon" }
	return "SiO2"
}

proc myName { x y z } { return 4e20 }

set size 50

param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I

init minx=0 miny=0 minz=0 maxx=24 maxy=$size maxz=$size material=material

profile name=I                proc=myName
profile name=C2 defect=CCluster proc=myName
insert interface particle=B coord={ 12 25 25 }

report defects

restart save=dump
