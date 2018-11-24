########## TEST 1
set sizeBox 20 
set temp -100 
set events 20000000 

set kB  8.6174e-5 
set l [expr [param get type=float key=S_Iron/Models/lambda]*1e-7] 
set TKelv [expr $temp+273.15]

proc material { x y z } { return "Iron" } 

param set type=map<string,string> key=MC/General/materials value={ Iron Fe }
param set type=bool   key=MC/Mesh/periodic.x value=true 

param set type=array<string,string> key=Iron/Models/interactions index=<111>+Gas value=false
param set type=array<string,string> key=Iron/Models/interactions index=VCluster+Gas value=false
param set type=map<string,bool>     key=Iron/Models/defined      index=HeCluster    value=true
param set type=map<string,bool>     key=Iron/Models/particles    index=He           value=true
param set type=map<string,bool>     key=Iron/Models/particles    index=HeV          value=true

param set type=array<string,string> key=Iron/Models/interactions value=false index=He+He 
param set type=array<string,string> key=Iron/Models/interactions value=false index=He+V 
param set type=array<string,string> key=Iron/Models/interactions value=false index=V+V 
param set type=array<string,string> key=Iron/Models/interactions value=false index=HeV+V

param set type=proc   key=Iron/HeCluster/migration value={ { return "HeV2 < 0% 1 0.09 20% 2 .07 100% 2 0.07 > " } } 

set PmHeV2 1
set EmHeV2 0.09

set CrPmHeV2 2
set CrEmHeV2 0.07

set DHeV2   [expr $PmHeV2  *exp(-$EmHeV2/($kB*$TKelv))]
set CrDHeV2 [expr $CrPmHeV2*exp(-$CrEmHeV2/($kB*$TKelv))]

#lowmsg "PmHeV2 is $PmHeV2, EmHeV2 is $EmHeV2, PbHeV2 is $PbHeV2 and EbHeV2 is $EbHeV2."

init material=material minx=0 miny=0 minz=0 maxx=$sizeBox maxy=$sizeBox maxz=$sizeBox 

for {set i 0} {$i<500} {incr i +1} {
	set ind [expr 1234+$i]	
	expr srand($ind)
	set x1 [expr rand()*$sizeBox]
	set y1 [expr rand()*$sizeBox]
	set z1 [expr rand()*$sizeBox]

	insert defect=HeCluster ID=HeV2 coord={$x1 $y1 $z1} no.print
}

anneal time=1e30 temp=$temp events=$events 

set r [extract diffusivity name=He material=Iron macroscopic] 
test tag=diffmic.HeV2 float=$r value=$DHeV2 error=0.05 

set density [param get key=Iron/Models/atomic.density type=float]

proc concent { x y z } {
	global density
	return [expr $density*0.45]
}


init material=material minx=0 miny=0 minz=0 maxx=$sizeBox maxy=$sizeBox maxz=$sizeBox 


for {set i 0} {$i<500} {incr i +1} {
	set ind [expr 1234+$i]	
	expr srand($ind)
	set x1 [expr rand()*$sizeBox]
	set y1 [expr rand()*$sizeBox]
	set z1 [expr rand()*$sizeBox]

	insert defect=HeCluster ID=HeV2 coord={$x1 $y1 $z1} no.print
}

profile proc=concent name=Cr

anneal time=1e30 temp=$temp events=$events 

set r [extract diffusivity name=He material=Iron macroscopic] 
test tag=diffmic.HeV2 float=$r value=$CrDHeV2 error=0.05 

