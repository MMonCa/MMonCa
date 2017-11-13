param set type=map<string,string>   key=MC/General/materials value="Gas Gas S_Iron Fe" 
set time 1
set conc 4e20

proc material { x y z } { 
	if { $x < 0 } { return "Gas" }
	return "S_Iron" 
}

proc myName { x y z } {
  global conc
  if { $x > 9 && $x < 12 } { return $conc }
  return 0
}

set size 30

param set type=map<string,float> key=Gas_S_Iron/Vacancy/barrier.right value=10 index=V
param set type=map<string,float> key=Gas_S_Iron/Iron/barrier.right    value=10 index=I
param set type=arrhenius         key=S_Iron/Vacancy/V(formation)    value={ 1e-5 10 }
param set type=arrhenius         key=S_Iron/Helium/HeV(formation)     value={ 3.93e-26 4.1 }
param set type=array<string,string> key=S_Iron/Models/interactions    value=false index=He+V
param set type=array<string,string> key=S_Iron/Models/interactions    value=false index=HeV+HeV
param set type=map<string,float>  key=Gas_S_Iron/Helium/desorption.low       value=0     index=HeV
param set type=map<string,float>  key=Gas_S_Iron/Helium/desorption.high      value=1     index=HeV
param set type=map<string,float>  key=Gas_S_Iron/Helium/desorption.threshold value=1e14  index=He
init minx=-1.5 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

profile name=HeV proc=myName

set howMany [expr $size*$size*3*1e-21 * $conc]
set trapped [expr $size*$size*1e-14*1e14]

test tag=before float=[extract count.particles particle=HeV] value=$howMany  error=0.05

#proc snapshot { } {
#	lowmsg "Counting [extract time] [extract count.particles particle=He defect=Interface] [extract count.particles]"
#}

anneal time=1e10 temp=700 events=3e7
#snapshot

test tag=after float=[extract count.particles] value=$trapped  error=.05
