param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"

set size 15
set time 300

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false
param set type=float  key=MC/Mesh/spacing.x value=0.5
param set type=float  key=MC/Mesh/spacing.y value=1.0
param set type=float  key=MC/Mesh/spacing.z value=1.0
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster:*+V       value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster:*+HeV     value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster:*+He      value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster:*+I       value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster+HeCluster value=true

param set type=array<string,string> key=S_Iron/Models/interactions index=He+He    value=HeCluster,1
param set type=array<string,string> key=S_Iron/Models/interactions index=He+HeV   value=HeCluster,1
param set type=array<string,string> key=S_Iron/Models/interactions index=HeV+HeV  value=HeCluster,1
param set type=array<string,string> key=S_Iron/Models/interactions index=He+V  value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=He+I  value=false

param set type=arrhenius key=S_Iron/Helium/He(migration) value={ 1 .4 }


param set type=proc   key=S_Iron/<111>/migration             value={ {
        set list ""
        for { set size 2 } { $size < 500 } { incr size } {
                set pref [expr 3.5e-4 + 1.7e-3/pow($size,1.7)]
                lappend list I$size < $pref 1.0 >
        }
	return $list
} }
param set type=string key=S_Iron/HeCluster/expand.impurity value=He new
param set type=float  key=S_Iron/HeCluster/expand.impurity.volume.nm3 value=0.26 new
param set type=float  key=S_Iron/HeCluster/expand.impurity.radius.M value=0.2 new
param set type=float  key=S_Iron/HeCluster/expand.impurity.radius.P value=0.5 new


init minx=-20 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material


proc snapshot { } { }

cascade file=cascades periodic fluence=3e14 flux=1e12 temp=-75 correct.for.surface

     set average 0
     set howmany 0
     foreach { i j k } [extract fuzz] {
        incr howmany
	set average [expr $average + $k]
	lowmsg "($i, $j) -> $k"
     }
     save lammps=nodist-evolution

     test tag=average float=-1.33 value=[expr $average/$howmany] error=.03
