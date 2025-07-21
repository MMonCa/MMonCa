param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas" 

param set type=float key=MC/Mesh/spacing.x value=6
param set type=float key=MC/Mesh/spacing.y value=6
param set type=float key=MC/Mesh/spacing.z value=6

set sizeX   60
set sizeY  300
set sizeZ  300
set T      900
set time    60
set Bcon  1e20

set kB  8.6174e-5

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}



proc spike { x y z } {
  global Bcon
  if { $x > 15 && $x < 45 } { return $Bcon }
  return 0
}

param set type=int key=MC/General/domains value=1

param set type=map<string,bool> key=Silicon/Models/particles index=I_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=I       value=true
param set type=map<string,int>  key=Silicon/Silicon/I(state.charge)        value="I 0"
param set type=map<string,bool> key=Silicon/Models/particles index=V_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=V       value=true
param set type=map<string,int>  key=Silicon/Vacancy/V(state.charge)        value="V 0"
param set type=map<string,bool> key=Silicon/Models/particles index=BI_0,-    value=false
param set type=map<string,bool> key=Silicon/Models/particles index=BI        value=true
param set type=map<string,int>  key=Silicon/Boron/BI(state.charge)           value="BI -1"

param set type=array<string,string> key=Silicon/Models/interactions index=B+BI  value=false
param set type=array<string,string> key=Silicon/Models/interactions index=I+I value=false

param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 9.21e29 3.74 }
param set type=arrhenius key=Silicon/Silicon/I(formation) value={ 6e26 4 }
param set type=arrhenius key=Silicon/Vacancy/V(migration) value={ 5.00e-8 0.4 } new
param set type=arrhenius key=Silicon/Silicon/I(migration) value={ 5e-2 0.8 }    new
param set type=arrhenius key=Silicon/Boron/BI(formation)    value={ 920 2.9 }
param set type=arrhenius key=Silicon/Boron/BI(migration)    value={ 5e-3 0.77 } new

set PfI   [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 0]
set EfI   [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set PfB   [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 0]
set EfB   [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 1]

set Icon [expr $PfI * exp (- $EfI / ($kB * ($T + 273.15)))]

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material

profile proc=spike name=B
profile proc=spike name=I
profile proc=spike name=V


anneal temp=$T time=$time

# gatherSpecificAtomFromInterfaceAndCluster
test array="[extract profile defect=* name=B]" value=[expr exp(log(10) * 19.1)] error=0.2 init=20 end=45 tag=allBclustered

# gatherSpecificAtomFromSpecificClusterFamily
test array="[extract profile defect=BICs material=Silicon name=B*]" value=[expr exp(log(10) * 19.1)] error=0.2 init=20 end=45 tag=allBinBICs

# gatherClusterFamily
test array="[extract profile defect=BICs material=Silicon]" value=[expr exp(log(10) * 18.75)] error=0.2 init=20 end=45 tag=allBICpieces

# gatherSpecificClusterFromClusterFamily
test array="[extract profile defect=BICs material=Silicon name=B3]" value=[expr exp(log(10) * 18.25)] error=0.25 init=20 end=45 tag=B3piecesInBIC

# gatherVanillaParticle
test array="[extract profile name=B]" value=[expr exp(log(10) * 20)] error=0.2 init=20 end=45 tag=allBvanilla

# gatherAllInActive
test array="[extract profile name=B*]" value=[expr exp(log(10) * 19.95)] error=0.2 init=20 end=45 tag=allBactive

