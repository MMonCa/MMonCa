param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

########## TEST 1
set sizeBox 20 
set temp -100 
set events 20000000 

set kB  8.6174e-5 
set l [expr [param get type=float key=S_Iron/Models/lambda]*1e-7] 
set TKelv [expr $temp+273.15]

proc material { x y z } { return "S_Iron" } 

param set type=bool   key=MC/Mesh/periodic.x value=true 
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=He+He 
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=He+V 
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=V+V 
param set type=array<string,string> key=S_Iron/Models/interactions value=HeCluster,1  index=HeV+V

param set type=arrhenius key=S_Iron/Helium/HeV(formation) value="1 1.008" 
param set type=arrhenius key=S_Iron/Vacancy/V(migration) value="1 0.25" 
param set type=arrhenius key=S_Iron/Helium/HeV(migration) value="0 0.25" 

set PmV    [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 0]
set EmV    [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 1]
set PbHeV  [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)]  0]
set EfHeV  [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)]  1]
set EfHe   [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 1]
set EfV    [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 2]
set EbHeV  [expr $EfHe + $EfV - $EfHeV]

param set type=proc   key=S_Iron/HeCluster/formation value={ { global EfV; return "HeV2 < [expr 2*$EfV-1.5] > " } }
param set type=proc   key=S_Iron/HeCluster/prefactor value={ { return "HeV2,V 1" } }
param set type=proc   key=S_Iron/HeCluster/migration value={ { return "HeV2 < 1 0.09 > " } } 

set PmHeV2 1
set EmHeV2 0.09

set DHeV2  [expr $PmHeV2*exp(-$EmHeV2/($kB*$TKelv))]

#lowmsg "PmHeV2 is $PmHeV2, EmHeV2 is $EmHeV2, PbHeV2 is $PbHeV2 and EbHeV2 is $EbHeV2."

init material=material minx=0 miny=0 minz=0 maxx=$sizeBox maxy=$sizeBox maxz=$sizeBox 

for {set i 0} {$i<500} {incr i +1} {
	set ind [expr 1234+$i]	
	expr srand($ind)
	set x1 [expr rand()*$sizeBox]
	set y1 [expr rand()*$sizeBox]
	set z1 [expr rand()*$sizeBox]

	insert defect=HeCluster ID=HeV2 coord={$x1 $y1 $z1}
}

anneal time=1e30 temp=$temp events=$events 

set r [extract diffusivity name=He material=S_Iron macroscopic] 
test tag=diffmic.HeV2 float=$r value=$DHeV2 error=0.05 

#### test 2

param set type=arrhenius key=S_Iron/Vacancy/V(migration) value="0.1 0.15" 
param set type=arrhenius key=S_Iron/Helium/HeV(migration)  value="0 0.45" 

param set type=proc   key=S_Iron/HeCluster/formation    value={ { global EfV; return "HeV2 < [expr 2*$EfV -1.5] > " } }
param set type=proc   key=S_Iron/HeCluster/prefactor    value={ { return "HeV2,V 5e1" } }
param set type=proc   key=S_Iron/HeCluster/migration    value={ { return "HeV2 < 1e-3 0.15 > " } } 

set fact     [expr 6./pow($l,2)]
set PmV      [expr [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 0]*$fact]
set EmV      [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 1]
set PmHeV2   [expr 1e-3*$fact]
set EmHeV2   0.15
set PHeV2    [expr 5e1*$fact]
set EpotHeV2 -1.5

#lowmsg "PmHeV2 is $PmHeV2, EmHeV2 is $EmHeV2, PHeV2 is $PHeV2 and EpotHeV2 is $EpotHeV2."

init material=material minx=0 miny=0 minz=0 maxx=$sizeBox maxy=$sizeBox maxz=$sizeBox 


for {set i 0} {$i<300} {incr i +1} {
	set ind [expr 1234+$i]	
	expr srand($ind)
	set x1 [expr rand()*$sizeBox]
	set y1 [expr rand()*$sizeBox]
	set z1 [expr rand()*$sizeBox]

	insert particle=HeV coord={$x1 $y1 $z1}
}

for {set j 0} {$j<300} {incr j +1} {
	set jnd [expr 5678+$j]	
	expr srand($jnd)
	set x2 [expr rand()*$sizeBox]
	set y2 [expr rand()*$sizeBox]
	set z2 [expr rand()*$sizeBox]

	insert particle=V coord={$x2 $y2 $z2}
}

set numHeVi   [extract count.particles particle=HeV]
set numVi     [extract count.particles particle=V]
set numHeV2i  [extract count.defects defect=Cluster ID=HeV2]
#lowmsg "#HeV es $numHeVi, #V es $numVi y #HeV2 es $numHeV2i."

set HeVconc  [expr $numHeVi/pow($sizeBox,3)*1e21]
set Vconc    [expr $numVi/pow($sizeBox,3)*1e21]
set HeV2conc [expr $numHeV2i/pow($sizeBox,3)*1e21]
set conc     [expr ($numHeVi+$numHeV2i)/pow($sizeBox,3)*1e21]
lowmsg "HeVconc es $HeVconc, Vconc es $Vconc, HeV2conc es $HeV2conc y conc es $conc."

anneal time=1e30 temp=$temp events=[expr $events*3] 

#lowmsg "EbHeV $EbHeV Epot $EpotHeV2 EmV $EmV"

set Ebarrier [expr -$EbHeV-$EpotHeV2+$EmV]
set break    [expr $PHeV2*exp(-$Ebarrier/($kB*$TKelv))]
set numHeV   [extract count.particles particle=HeV]
set numV     [extract count.particles particle=V]
set numHeV2  [extract count.defects defect=Cluster ID=HeV2]
set vol      [expr 4./3.*3.1415926*pow($l,3)]
set migV     [expr $PmV*exp(-$EmV/($kB*$TKelv))]
set C        [expr $break/($vol*$migV)]

#lowmsg "Ebarrier es $Ebarrier, PHeV2 es $PHeV2, break es $break, vol es $vol and migV es $migV. Lambda es $l."
#lowmsg "el radicando vale: [expr pow($C,2)] mas [expr 4*$conc*$C] y es igual a [expr (pow($C,2)+4*$conc*$C)]"
#lowmsg "La raiz es [expr sqrt(pow($C,2)+4*$conc*$C)]"
#lowmsg "C es $C."

set rest [expr (1./2.*($C-sqrt($C*($C+4*$conc))))]
set HeV2conc  [expr ($conc+$rest)]

#lowmsg "Rest es $rest y HeV2conc es $HeV2conc"

set Ratio  [expr $HeV2conc/$conc]
set DHeV2  [expr $PmHeV2*exp(-$EmHeV2/($kB*$TKelv))/$fact]
set DHe    [expr $Ratio*$DHeV2]

#lowmsg "Ratio is $Ratio"
#lowmsg "DHeV2 es $DHeV2 y DHe es $DHe."

set r [extract diffusivity macroscopic name=He material=S_Iron]

test tag=diffmac.HeV2 float=$r value=$DHe error=0.05
