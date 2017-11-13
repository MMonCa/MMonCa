param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas" 

set time 1
set EV0 [expr 9998583*1.2]
set EV1 [expr 1417*1.2]

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

set size 20

param set type=bool      key=MC/Mesh/speed.up  value=false
param set type=arrhenius key=S_Iron/Iron/I(formation)   value={ 1e-5 10 }
param set type=arrhenius key=S_Iron/Iron/I(migration)   value={ 1.0e-3 1.0 }
param set type=arrhenius key=S_Iron/Vacancy/V(formation)        value={ 1e-5 10 }
param set type=array<string,string> key=S_Iron/Models/interactions index=I+Gas value=false
param set type=array<string,string> key=S_Iron/Models/interactions index=V+Gas value=false
param set type=proc      key=S_Iron/ICluster/formation  value={ {  
	set list ""
	set eF 10
	for { set size 2 } { $size < 9 } { incr size } {
		if { $size < 4 } {  
			lappend list I$size
			lappend list < [expr $size*$eF + 0] >
		} else {
			lappend list I$size
			lappend list < [expr $size*$eF + 0.4] >
		}
       }
       return $list
} }
param set type=proc      key=S_Iron/ICluster/prefactor  value={ {                                  
        return "I4,I 1e-2"
} }
param set type=proc      key=S_Iron/ICluster/migration  value={ {
        set list ""
        for { set size 2 } { $size < 9 } { incr size } {
             lappend list I$size < 0 5 >
       }
       return $list
} }


foreach { T } { 1100 900 700 500 300 } {

init minx=-2 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=ICluster ID=I4 coord={ [expr $size/2] [expr $size/2] [expr $size/2] }

set kB  8.6174e-5
set l         [param get type=float               key=S_Iron/Models/lambda]
set I_P       [lindex [param get type=arrhenius   key=S_Iron/Iron/I(migration)] 0]
set I_E       [lindex [param get type=arrhenius   key=S_Iron/Iron/I(migration)] 1]
set Cluster_P 1e-2
set Cluster_E 0.4

set convFactor [expr 6./$l/$l*1e14]
set time0 [expr 1./($I_P*$convFactor*exp(-$I_E/($kB*($T+273.15))))]
set time1 [expr 1./($Cluster_P*$convFactor*exp(-($Cluster_E+$I_E)/($kB*($T+273.15))))]
set time  [expr $EV0*$time0 + $EV1*$time1]
set CI    [expr $EV0/$size/$size/$size/$time*1e21*$time0]

anneal time=1e25 temp=$T events=[expr $EV0+$EV1]

lowmsg "Time 0 $time0"
lowmsg "Time 1 $time1"

lowmsg "time should be $time"
lowmsg "Conc should be $CI"
lowmsg "be aware of the Moire effects..."
set average [test tag=$T.I array={ [extract profile.mobile name=I] } value=$CI init=0 end=$size error=.15]
test tag=$T.I.average float=$average value=$CI error=.05
test tag=$T.I.time    float=[extract time] value=$time error=.05
}
