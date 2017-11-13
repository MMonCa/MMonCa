param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas" 
set time 1
set conc 4e20

proc material { x y z } { 
	if { $x < 0 } { return "Gas" }
	return "S_Iron" 
}

set size 30

param set type=map<string,float>    key=Gas_S_Iron/Vacancy/barrier.right value=10 index=V
param set type=map<string,float>    key=Gas_S_Iron/Iron/barrier.right    value=10 index=I
param set type=arrhenius            key=S_Iron/Vacancy/V(formation)    value={ 1e-5 10 }
param set type=arrhenius            key=S_Iron/Iron/I(formation)    value={ 1e-5 10 }
param set type=arrhenius            key=S_Iron/Helium/HeV(formation)     value={ 3.93e-26 4.1 }
param set type=array<string,string> key=S_Iron/Models/interactions    value=false index=He+V
param set type=array<string,string> key=S_Iron/Models/interactions    value=false index=He+I
param set type=array<string,string> key=S_Iron/Models/interactions    value=false index=HeV+HeV
param set type=map<string,float>    key=Gas_S_Iron/Helium/desorption.low       value=0     index=HeV
param set type=map<string,float>    key=Gas_S_Iron/Helium/desorption.high      value=0     index=HeV
param set type=map<string,float>    key=Gas_S_Iron/Helium/desorption.threshold value=1e20  index=He
param set type=arrhenius            key=Gas_S_Iron/Helium/He(migration) value={ 1e-2 .5 }

init minx=-1.5 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert particle=HeV coord={ 10 15 15 }

test tag=before float=1 value=[extract count.particles particle=HeV]  error=0

anneal time=1e-5 temp=700

lowmsg [extract defects={ }]

report defects

test tag=after float=[extract count.particles particle=He] value=1  error=0

extract defects={ } file=nodist-before.txt

anneal time=1e-5 temp=700

extract defects={ } file=nodist-after.txt

set rc [catch {exec diff nodist-before.txt nodist-after.txt} output]
if {$rc == 0 } {
        test tag=comparison one=passed not equal two=passed
} else {
         test tag=comparison one=passed equal two=passed
}
