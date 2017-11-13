param set type=map<string,string>   key=MC/General/materials value="Gas Gas Silicon Si"

set size 10
set time 300

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false

param set type=float  key=Silicon/Models/amorphization.threshold  value=-1

param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+V

init minx=0 miny=0 minz=0 maxx=25 maxy=$size maxz=$size material=material

for { set i 0 } { $i < 5000 } { incr i } {


# put 10 cascades that is 5*5 + 5*7 atoms = 60
cascade file=data/cascades periodic fluence=1e13 temp=-250
lowmsg "Cascade $i"
}
report all

test tag=count   float=[extract count.particles]            value=300000 error=0
test tag=count.I float=[extract count.particles particle=I] value=175000 error=0
test tag=count.V float=[extract count.particles particle=V] value=125000 error=0
