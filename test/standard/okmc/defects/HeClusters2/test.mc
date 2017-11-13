param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 

#This test works on the assumption that the time for breaking H2V2 is going to be so high that it does not
#matter how much is the time for V or HeV diffusion
#But if the HeV2 prefactor is small, then the assumption is not correct

set time 1
set ev_He(1000)    14271731
set ev_V2(1000)    36000000
set ev_HeV2(1000)  1367

proc material { x y z } { return "S_Iron" }

set size 20

param set type=bool   		    key=MC/Mesh/speed.up  	    value=false
param set type=arrhenius            key=S_Iron/Helium/He(migration) value="1e-3 1.0"
param set type=array<string,string> key=S_Iron/Models/interactions  value=true    index=HeCluster:V2+He 
param set type=array<string,string> key=S_Iron/Models/interactions  value=false   index=V+V
param set type=array<string,string> key=S_Iron/Models/interactions  value=false   index=V+He
param set type=array<string,string> key=S_Iron/Models/interactions  value=false   index=I+He

set kB  8.6174e-5
set l         [param get type=float             key=S_Iron/Models/lambda]
set He_Pf     [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 0]
set He_Ef     [lindex [param get type=arrhenius key=S_Iron/Helium/He(formation)] 1]
set He_Pm     [lindex [param get type=arrhenius key=S_Iron/Helium/He(migration)] 0]
set He_Em     [lindex [param get type=arrhenius key=S_Iron/Helium/He(migration)] 1]
set V2_Pm    1e-3
set V2_Em    0.9
set HeV2_Ef  9
set V2_Ef    9.5

set capVol    [expr 4./3.*3.1415926*$l*$l*$l*0.98*0.98*0.98*1e-21]
set convFactor [expr 6./$l/$l*1e14]

param set type=proc key=S_Iron/HeCluster/prefactor value={ {
        set list ""
        append list     "HeV2,He 1e-18 "
        return $list
} }

param set type=proc key=S_Iron/HeCluster/migration value={ {
        set list ""
        append list     "V2 < 1e-3 0.9 >"
        return $list
} }

param set type=proc key=S_Iron/HeCluster/formation value={ {
        set list ""
        append list "V2   < 9.5 > "
        append list "HeV2 < 9 > "
        return $list
} }

foreach { T } { 1000 } {

set EV_He   $ev_He($T)
set EV_V2   $ev_V2($T)
set EV_HeV2 $ev_HeV2($T)

set freqHem [expr $He_Pm*$convFactor*exp(-$He_Em/($kB*($T+273.15)))]

set Cluster_P   [expr $convFactor*1e-18]
set energ       [expr ($He_Ef + $V2_Ef - $HeV2_Ef >= 0)? $He_Ef + $V2_Ef - $HeV2_Ef : 0]
set freqHeV2_V2 [expr $Cluster_P*exp(-($energ+$V2_Em)/($kB*($T+273.15)))]
set freqHeV2_He [expr $Cluster_P*exp(-($energ+$He_Em)/($kB*($T+273.15)))]
set timeHeV2    [expr 1./($freqHeV2_He + $freqHeV2_V2)]

set time [expr $EV_HeV2*$timeHeV2]
set CHe  [expr $EV_He  /$size/$size/$size/$time*1e21/$freqHem]

lowmsg "time should be $time"
lowmsg "CHe  should be $CHe"

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
insert defect=HeCluster ID=HeV2 coord={ [expr $size/2] [expr $size/2] [expr $size/2] }
anneal time=1e25 temp=$T events=5e7

lowmsg "be aware of the Moire effects..."
test              tag=$T            float=[extract time]                       value=$time                 error=.1
set average [test tag=$T.He         array={ [extract profile.mobile name=He]} value=$CHe init=0 end=$size error=.15]
test              tag=$T.He.average float=$average                             value=$CHe                  error=.1
}
