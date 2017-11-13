param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

set size 20
set temp -100
set events 15000000
set conc 1e20

set kB  8.6174e-5
proc constant { x y z } { global conc; return $conc }
proc material { x y z } { return "S_Iron" }
param set type=bool   key=MC/Mesh/periodic.x value=true

param set type=array<string,string> key=S_Iron/Models/interactions value=false index=He+He
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=HeV+V

#---------------------------------------------------

param set type=array<string,string> key=S_Iron/Models/interactions value=true index=He+V
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=V+V

param set type=arrhenius key=S_Iron/Helium/He(migration) value="0 5"
param set type=arrhenius key=S_Iron/Helium/HeV(migration)  value="1 .25"
param set type=arrhenius key=S_Iron/Vacancy/V(migration) value="1 .25"
param set type=arrhenius key=S_Iron/Vacancy/V(formation) value="5e-3 1.258"
param set type=arrhenius key=S_Iron/Helium/HeV(formation)  value="1 -10"

set l [expr [param get type=float key=S_Iron/Models/lambda]*1e-7]
set PmHeV [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(migration)]  0]
set EmHeV [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(migration)]  1]
set PmV   [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 0]
set EmV   [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 1]

set DHeV  [expr $PmHeV*exp(-$EmHeV/($kB*($temp+273.15)))]
set v_cap [expr 3.65*$l*$l*$l]
set form  [expr 4./3.*3.1415926*$l*$l*$l *$PmV* exp(-$EmHeV/($kB*($temp+273.15)))]

#-------- 

init material=material minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size

profile proc=constant name=HeV

anneal time=1e30 temp=$temp events=$events

set r [extract diffusivity name=HeV material=S_Iron]
test tag=diff.HeV float=$r value=$DHeV error=0.05

#-------
param set type=arrhenius key=S_Iron/Helium/HeV(formation) value="5.85e-27 1.208"
set PfHeV [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)] 0]
set EfHeV [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)] 1]
set PfV   [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 0]
set EfV   [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 1]
set PfHe  [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 0]
set EfHe  [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 1]
set PbHeV [expr $PfHe * $PfV * $PmV * $v_cap / $PfHeV]
set EbHeV [expr $EfHe + $EfV - $EfHeV]
set break [expr $PbHeV * exp(-($EmV+$EbHeV)/($kB*($temp+273.15)))]

init material=material minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size

profile proc=constant name=He
profile proc=constant name=V

anneal time=1e30 temp=$temp events=$events

lowmsg " Vconc [extract count.particles particle=V].0/($size*$size*$size)*1e21"
set Vconc [expr [extract count.particles particle=V].0/($size*$size*$size)*1e21]
set r [extract diffusivity macroscopic name=He material=S_Iron]
lowmsg "ratio $Vconc*$form/($Vconc*$form+$break)"
set ratio [expr $Vconc*$form/($Vconc*$form+$break)]
set DHe [expr $ratio*$DHeV]

lowmsg "measured $r, Theoretical $DHe, from $ratio and $DHeV diffusivity, $Vconc is V"
test tag=diff.brHeV float=$r value=$DHe error=0.08

#------
param set type=arrhenius key=S_Iron/Helium/He(migration) value="1e-3 .5"

init material=material minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size

profile proc=constant name=He

anneal time=1e30 temp=$temp events=$events

set r [extract diffusivity macroscopic name=He material=S_Iron]
set PmHe [lindex [param get type=arrhenius key=S_Iron/Helium/He(migration)] 0]
set EmHe [lindex [param get type=arrhenius key=S_Iron/Helium/He(migration)] 1]
set D [expr $PmHe*exp(-$EmHe/($kB*($temp+273.15)))]
test tag=diff.He float=$r value=$D error=0.08
