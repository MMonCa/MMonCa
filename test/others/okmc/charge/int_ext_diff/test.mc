set size 20
set T 600
set events 5e7
set conc 6e18

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"

# activar carga y ver que pasa!

param set type=bool  key=MC/Electrostatic/load.poisson value=true
param set type=float key=MC/Electrostatic/update.events value=1e6
param set type=int   key=MC/Electrostatic/update.time.decade value=3

#param set type=map<string,bool> key=Silicon/Models/particles index=Bi_0,-  value=false
#param set type=map<string,bool> key=Silicon/Models/particles index=I_0,-,+ value=false
#param set type=map<string,bool> key=Silicon/Models/particles index=Bi_-  value=true
#param set type=map<string,bool> key=Silicon/Models/particles index=I_0   value=true
param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 1e5 8 }

set kB  8.6174e-5
proc constant  { x y z } { global conc; return $conc }
proc constant2 { x y z } { global conc; return [expr $conc/100] }
proc material { x y z } { if {$x<0} { return "Gas" }; return "Silicon" }

param set type=bool   key=MC/Mesh/periodic.x value=false

param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+BI
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+BI
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+I
param set type=array<string,string> key=Silicon/Models/interactions value=true  index=B+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I

set If_P [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set If_E [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 2]
set CI0  [expr $If_P*exp(-$If_E/($kB*($T+273.15)))]

set Bif_P [lindex [param get type=arrhenius key=Silicon/Boron/BI(formation)] 0]
set Bif_E [lindex [param get type=arrhenius key=Silicon/Boron/BI(formation)] 1]

set BimM_P [lindex [param get type=arrhenius key=Silicon/Boron/BI_-(migration)] 0]
set BimM_E [lindex [param get type=arrhenius key=Silicon/Boron/BI_-(migration)] 1]
set Bim0_P [lindex [param get type=arrhenius key=Silicon/Boron/BI_0(migration)] 0]
set Bim0_E [lindex [param get type=arrhenius key=Silicon/Boron/BI_0(migration)] 1]

set eBiM    [param get type=float key=Silicon/Boron/BI(e(-1,0))]

set Bf_P  [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 1]  
set Bf_E  [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 2]  

foreach newconc { 1e19 5e18 2e19 } {

set conc $newconc

init material=material minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size

profile proc=constant name=B
set CB [extract profile name=B dimension=0]
lowmsg "Initial B conc is $CB versus $conc" 

#profile proc=constant2 name=I

anneal time=1e30 temp=$T events=$events

set CI   [extract profile.mobile name=I state=0 dimension=0]
test tag=$T.$conc.I float=$CI value=$CI0 error=.25

#set CBi    [expr $CI/$CI0*$CB*$Bif_P*exp(-$Bif_E/($kB*($T+273.15)))]

set Eg      [extract electrostatic name=bandgap dimension=0 min.x=0 max.x=$size]
set Eg0     [Silicon/ElectronicStructure/Bandgap 300]
set scale   [expr $Eg/$Eg0]
lowmsg      "Eg is $Eg, initial is $Eg0, scale is $scale"
set CBiM    [expr $CB*$Bif_P*exp(-$Bif_E/($kB*($T+273.15)))]
set eF      [expr -[extract electrostatic name=valence.bandedge dimension=0 min.x=0 max.x=$size] ]
lowmsg      "eF is $eF and eBiM is $eBiM"
set CBi0    [expr $CBiM*exp(($scale*$eBiM - $eF)/($kB*($T+273.15)))]

test tag=$T.$conc.hole array={[extract electrostatic  name=hole.density dimension=1]} value=$conc error=.76 init=5 end=$size
test tag=$T.$conc.CBiM array={[extract profile.mobile name=BI state=- dimension=1]} value=$CBiM error=.60 init=5 end=$size
test tag=$T.$conc.CBi0 array={[extract profile.mobile name=BI state=0 dimension=1]} value=$CBi0   error=.75 init=5 end=$size

set DBM     [expr $BimM_P * exp(-$BimM_E/($kB*($T+273.15)))]
set DB0     [expr $Bim0_P * exp(-$Bim0_E/($kB*($T+273.15)))]
set thDiff  [expr $DBM*$CBiM/$CB + $DB0*$CBi0/$CB]

test tag=$T.$conc.diff float=[extract diffusivity macroscopic material=Silicon name=B] value=$thDiff error=.75
}
