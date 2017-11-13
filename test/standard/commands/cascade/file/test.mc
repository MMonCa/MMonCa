param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=float		    key=Silicon/Models/amorphization.threshold value=-1
param set type=bool                 key=MC/General/average.rates value=true

set size 10
set time 300

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false

param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+V

init minx=0 miny=0 minz=0 maxx=25 maxy=$size maxz=$size material=material

#put 10 cascades that is 5*5 + 5*7 atoms = 60
cascade file=data/cascades periodic fluence=1e13 flux=1e12 temp=-250

test tag=count   float=[extract count.particles]              value=60 error=0.1
test tag=count.I float=[extract count.particles particle=I] value=35 error=0.1
test tag=count.V float=[extract count.particles particle=V] value=25 error=0.1

report all

extract defects={ } file=p.txt
exec cut -d " " -f 1-6 <p.txt >nodist-particles.txt
exec rm p.txt
exec diff data/particles.txt nodist-particles.txt

set size 100
init minx=0 miny=0 minz=0 maxx=25 maxy=$size maxz=$size material=material

#put approx 10 cascades that is 500*5 + 500*7 atoms = 6000
cascade file=data/cascades periodic fluence=1e13 flux=1e12 temp=-250

report all

test tag=count   float=[extract count.particles]              value=6000 error=0.075
test tag=count.I float=[extract count.particles particle=I] value=3500 error=0.075
test tag=count.V float=[extract count.particles particle=V] value=2500 error=0.075
