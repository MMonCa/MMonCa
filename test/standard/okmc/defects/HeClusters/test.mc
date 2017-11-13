param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

#This test works on the assumption that the time for breaking H2V2 is going to be so high that it does not
#matter how much is the time for V or HeV diffusion
#But if the HeV2 prefactor is small, then the assumption is not correct

param set type=proc key=S_Iron/HeCluster/prefactor value={ {
	set list ""	
	append list     "He2V,HeV  1e-3 "
	append list	"He2V,V    1e-3 "
	append list	"HeV2,HeV  1e-18 "
	append list	"HeV2,V    1e-18 "
	append list	"He2V2,HeV 1e-3 "
	append list	"He2V2,V   1e-3"
	
	return $list
} }

param set type=proc key=S_Iron/HeCluster/formation value={ {
        set list ""
	set eF 10
 	append list "He2   < 2 > "
 	append list "He2V  < [expr 2*$eF -1] > "
 	append list "HeV2  < [expr 2*$eF -2] > "
 	append list "He2V2 < [expr 2*$eF -1.5] > "
        return $list
} }

set time 1
set ev_V(1000)    22267939
set ev_HeV(1000)  77252738
set ev_HeV2(1000) [expr 117+736]

proc material { x y z } { return "S_Iron" }

set size 20

param set type=bool   		     key=MC/Mesh/speed.up  	    value=false
param set type=arrhenius             key=S_Iron/Vacancy/V(formation)  value={ 1e-5 10 }
param set type=arrhenius             key=S_Iron/Vacancy/V(migration)  value={ 1.0e-3 1.0 }
param set type=arrhenius             key=S_Iron/Helium/HeV(formation) value={ 5.85e-30 9.1 }
param set type=arrhenius             key=S_Iron/Iron/I(formation)     value={ 1e-5 10 }
param set type=arrhenius             key=S_Iron/Helium/HeI(formation) value={ 2.4e-29 9.2 }
param set type=array<string,string>  key=S_Iron/Models/interactions   value=false   index=V+V
param set type=array<string,string>  key=S_Iron/Models/interactions   value=false   index=V+Gas
param set type=array<string,string>  key=S_Iron/Models/interactions   value=false   index=I+Gas

foreach { T } { 1000 } {

set EV_V    $ev_V($T)
set EV_HeV  $ev_HeV($T)
set EV_HeV2 $ev_HeV2($T)

set kB  8.6174e-5
set l         [param get type=float             key=S_Iron/Models/lambda]
set V_Pm      [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 0]
set V_Em      [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(migration)] 1]
set V_Pf      [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 0]
set V_Ef      [lindex [param get type=arrhenius key=S_Iron/Vacancy/V(formation)] 1]
set He_Pf     [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 0]
set He_Ef     [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 1]
set HeV_Pm    [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(migration)] 0]
set HeV_Em    [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(migration)] 1]
set HeV_Pf    [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)] 0]
set HeV_Ef    [lindex [param get type=arrhenius key=S_Iron/Helium/HeV(formation)] 1]

set capVol    [expr 4./3.*3.1415926*$l*$l*$l*0.98*0.98*0.98*1e-21]
set HeV_Eb    [expr $He_Ef + $V_Ef - $HeV_Ef]
set HeV_Pb    [expr $He_Pf * $V_Pf * $V_Pm * $capVol / $HeV_Pf]


set convFactor [expr 6./$l/$l*1e14]
set freqV     [expr $V_Pm*$convFactor*exp(-$V_Em/($kB*($T+273.15)))]
set freqHeVm  [expr $HeV_Pm*$convFactor*exp(-$HeV_Em/($kB*($T+273.15)))]
set freqHeVb  [expr $HeV_Pb*$convFactor*exp(-$HeV_Eb/($kB*($T+273.15)))]

set Cluster_P    [expr $convFactor*1e-18]
set Cluster_E    -2

set energ        [expr ($Cluster_E + $HeV_Eb >= 0)? 0 : -$Cluster_E - $HeV_Eb]
set freqHeV2_V   [expr $Cluster_P*exp(-($energ+$V_Em)  /($kB*($T+273.15)))]
set freqHeV2_HeV [expr $Cluster_P*exp(-($energ+$HeV_Em)/($kB*($T+273.15)))]
set timeHeV2     [expr 1./($freqHeV2_V + $freqHeV2_HeV)]

set time     [expr $EV_HeV2*$timeHeV2]
set CV       [expr $EV_V  /$size/$size/$size/$time*1e21/$freqV]
set CHeV     [expr $EV_HeV/$size/$size/$size/$time*1e21/$freqHeVm]

lowmsg "HeV binding is $HeV_Eb"
lowmsg "time should be $time"
lowmsg "CV   should be $CV"
lowmsg "CHeV should be $CHeV"

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
insert defect=HeCluster ID=HeV2 coord={ [expr $size/2] [expr $size/2] [expr $size/2] }
anneal time=1e25 temp=$T events=1e8

lowmsg "be aware of the Moire effects..."
test              tag=$T             float=[extract time]                       value=$time                  error=.12
set average [test tag=$T.V           array={ [extract profile.mobile name=V] }  value=$CV  init=0 end=$size  error=.15]
test              tag=$T.V.average   float=$average                             value=$CV                    error=.12
set average [test tag=$T.HeV         array={ [extract profile.mobile name=HeV]} value=$CHeV init=0 end=$size error=.15]
test              tag=$T.HeV.average float=$average                             value=$CHeV                  error=.12
}
