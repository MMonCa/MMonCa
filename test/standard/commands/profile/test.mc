param set type=map<string,string>   key=MC/General/materials value="Silicon Si SiO2 SiO2"
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

report all

test array="[extract profile name=I]" value=4e20 error=0.06  init=1 end=11 tag=I
test array="[extract profile name=C]" value=8e20 error=0.06  init=1 end=11 tag=C
test array="[extract profile name=I]" value=0 error=0.  init=13 end=24 tag=I
test array="[extract profile name=C]" value=0 error=0.  init=13 end=24 tag=C

test float=[extract count.defects ID=C2 defect=CCluster] value=[expr 12*$size*$size*1e-21*4e20] error=0.05 tag=Defects

param set type=int key=MC/General/domains value=2

init minx=0 miny=0 minz=0 maxx=24 maxy=$size maxz=$size material=material

profile name=I                proc=myName
profile name=C2 defect=CCluster proc=myName

test array="[extract profile name=I]" value=4e20 error=0.06  init=1 end=11 tag=I.par
test array="[extract profile name=C]" value=8e20 error=0.06  init=1 end=11 tag=C.par
test array="[extract profile name=I]" value=0 error=0.  init=13 end=24 tag=I.par
test array="[extract profile name=C]" value=0 error=0.  init=13 end=24 tag=C.par

test float=[extract count.defects ID=C2 defect=CCluster] value=[expr 12*$size*$size*1e-21*4e20] error=0.05 tag=Defects
