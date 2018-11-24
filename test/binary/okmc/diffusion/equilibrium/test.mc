set time 1

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "SiliconCarbide"
}

set sizeX  8
set sizeYZ 25

param set type=map<string,string> key=MC/General/materials value={ SiliconCarbide SiC Gas Gas }
param set type=array<string,string> key=SiliconCarbide/Models/interactions value={ {
	SiI+Gas true
	VSi+Gas true
	CI+Gas  true
	VC+Gas  true
} }

#param set type=map<string,bool> key=SiliconCarbide/Models/particles index=VC_0,-,-2,+,+2 value=false
#param set type=map<string,bool> key=SiliconCarbide/Models/particles index=VC_0 value=true

foreach { T } { 6000 3000 } {

init minx=-2 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

anneal time=1e25 temp=$T events=1e8

set kB  8.6174e-5

set VSiP [lindex [param get type=arrhenius key=SiliconCarbide/Vacancy/VSi(formation)] 0]
set VSiE [lindex [param get type=arrhenius key=SiliconCarbide/Vacancy/VSi(formation)] 1]
set VSiC [expr $VSiP*exp(-$VSiE/($kB*($T+273.15)))]
test tag=$T.VSi array={ [extract profile.mobile name=VSi] } value=$VSiC init=1 end=$sizeX error=0.15

set SiIP [lindex [param get type=arrhenius key=SiliconCarbide/Silicon/SiI(formation)] 0]
set SiIE [lindex [param get type=arrhenius key=SiliconCarbide/Silicon/SiI(formation)] 1]
set SiIC [expr $SiIP*exp(-$SiIE/($kB*($T+273.15)))]
#very noisy
#test tag=$T.SiI array={ [extract profile.mobile name=SiI] } value=$SiIC init=1 end=$sizeX error=0.15


set VCP [lindex [param get type=arrhenius key=SiliconCarbide/Vacancy/VC(formation)] 0]
set VCE [lindex [param get type=arrhenius key=SiliconCarbide/Vacancy/VC(formation)] 1]
set VCC [expr $VCP*exp(-$VCE/($kB*($T+273.15)))]
test tag=$T.VC array={ [extract profile.mobile name=VC] } value=$VCC init=1 end=$sizeX error=0.15

set CIP [lindex [param get type=arrhenius key=SiliconCarbide/Carbon/CI(formation)] 0]
set CIE [lindex [param get type=arrhenius key=SiliconCarbide/Carbon/CI(formation)] 1]
set CIC [expr $CIP*exp(-$CIE/($kB*($T+273.15)))]
test tag=$T.CI array={ [extract profile.mobile name=CI] } value=$CIC init=1 end=$sizeX error=0.15


}
