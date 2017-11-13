param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas" 
set time 1

proc material { x y z } { 
	if { $x < 0 } { return "Gas" }
	return "S_Iron" 
}

proc myName { x y z } {
  if { $x > 9 && $x < 12 } { return 2e20 }
  return 0
}

set size 30


param set type=map<string,float>     key=Gas_S_Iron/Vacancy/barrier.right  value=10        index=V
param set type=map<string,float>     key=Gas_S_Iron/Iron/barrier.right     value=10        index=I
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=He+V
param set type=array<string,string>  key=S_Iron/Models/interactions        value=false     index=HeV+HeV
param set type=array<string,string>  key=S_Iron/Models/interactions        value=true      index=HeCluster+Gas
param set type=proc                  key=S_Iron/HeCluster/migration        value={ { return "He2 < 2 0.8 > " } } 
param set type=map<string,float>     key=Gas_S_Iron/Helium/desorption.high value=1         index=HeCluster

init minx=-1.5 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

profile defect=HeCluster name=He2 proc=myName

set howMany [expr $size*$size*3*1e-21 * 2e20 * 2]
test tag=before float=[extract count.particles ID=He2] value=$howMany  error=0.05

anneal time=1e10 temp=700 events=1e7

test tag=after float=[extract count.particles] value=0  error=0

report reactions.interface
