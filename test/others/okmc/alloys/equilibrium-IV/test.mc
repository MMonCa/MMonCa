proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}
param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
set densityCm3 [param get type=float key=Silicon/Models/alloy.density]
set conc 0.2

proc prof { x y z } {
    global densityCm3
    global conc
    if { $x > 7 } {
    return [expr $conc * $densityCm3];
    }
    return 0
}

set sizeX  15
set sizeYZ 15
set T 1100

param set type=bool      key=MC/Mesh/speed.up value=false

param set type=map<string,bool> key=Silicon/Models/particles index=I_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=I       value=true
param set type=map<string,int>  key=Silicon/Silicon/I(state.charge)        value="I 0"
param set type=map<string,bool> key=Silicon/Models/particles index=V_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=V       value=true
param set type=map<string,int>  key=Silicon/Vacancy/V(state.charge)        value="V 0"

param set type=array<string,string> key=Silicon/Models/interactions index=I+I value=false
param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 1e-5 10 }
param set type=arrhenius key=Silicon/Vacancy/V(migration) value={ 1e-5 10 } new
param set type=arrhenius key=Silicon/Silicon/I(formation) value={ 0% 7.95e27 3.2 100% 4e26 4 }
param set type=arrhenius key=Silicon/Silicon/I(migration) value={ 0% 2e-3 0.7 100% 2e-3 0.791 } new

set IP_Si [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set IE_Si [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 2]
set IP_Ge [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 4]
set IE_Ge [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 5]

set IE_SiGe [expr (1 - $conc) * $IE_Si + $conc * $IE_Ge]
set IP_SiGe [expr ($IP_Si) ** (1 - $conc) * ($IP_Ge) ** $conc]

init minx=-1.5 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

profile name=Ge proc=prof

# Annealing for interstitials

anneal time=1e25 temp=$T events=2e7

set kB  8.6174e-5

set IC_Si [expr $IP_Si*exp(-$IE_Si/($kB*($T+273.15)))]
set IC_SiGe [expr $IP_SiGe*exp(-$IE_SiGe/($kB*($T+273.15)))]

lowmsg "expr $IP_SiGe*exp(-$IE_SiGe/($kB*($T+273.15))) is $IC_SiGe"

test tag=I.in.Si array={ [extract profile.mobile name=I file=I.profile] } value=$IC_Si init=0 end=7 error=0.05
test tag=I.in.SiGe array={ [extract profile.mobile name=I file=I.profile] } value=$IC_SiGe init=9 end=15 error=0.05

# Annealing for vacancies

param set type=array<string,string> key=Silicon/Models/interactions index=V+V value=false
param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 0% 7.95e27 3.2 100% 4e26 4 } 
param set type=arrhenius key=Silicon/Vacancy/V(migration) value={ 0% 2e-3 0.7 100% 2e-3 0.791 }
param set type=arrhenius key=Silicon/Silicon/I(formation) value={ 1e-5 10 }
param set type=arrhenius key=Silicon/Silicon/I(migration) value={ 1e-5 10 }

set VP_Si [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 1]
set VE_Si [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 2]
set VP_Ge [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 4]
set VE_Ge [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 5]

set VE_SiGe [expr (1 - $conc) * $VE_Si + $conc * $VE_Ge]
set VP_SiGe [expr ($VP_Si) ** (1 - $conc) * ($VP_Ge) ** $conc]

init minx=-1.5 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

profile name=Ge proc=prof

anneal time=1e25 temp=$T events=2e7

set kB  8.6174e-5

set VC_Si [expr $VP_Si*exp(-$VE_Si/($kB*($T+273.15)))]
set VC_SiGe [expr $VP_SiGe*exp(-$VE_SiGe/($kB*($T+273.15)))]

lowmsg "expr $VP_SiGe*exp(-$VE_SiGe/($kB*($T+273.15))) is $VC_SiGe"

test tag=V.in.Si array={ [extract profile.mobile name=V material=Silicon file=V.profile] } value=$VC_Si init=1 end=7 error=0.05
test tag=V.in.SiGe array={ [extract profile.mobile name=V material=Silicon] } value=$VC_SiGe init=9 end=15 error=0.05

