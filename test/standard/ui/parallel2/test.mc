param set type=map<string,string>   key=MC/General/materials value="S_Iron SFe Gas Gas"
param set type=int key=MC/General/domains    value=4
#param set type=int key=MC/General/subdomains value=2
#param set type=int key=MC/General/levels     value=3

set size 50
set time 300

proc material { x y z } { 
	if { $x < 0 } { return "Gas" }
	return "S_Iron" 
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false

for { set c 38 } { $c < 41 } { incr c } {

	expr srand($c)

	init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

	#introduce 100 FP and anneal a little...
	set time 1e-3
	set temp 75

	for { set i 0 } { $i < 100 } { incr i} {

	    insert particle=I coord={ [expr rand()*$size] [expr rand()*$size] [expr rand()*$size] }
	    insert particle=V coord={ [expr rand()*$size] [expr rand()*$size] [expr rand()*$size] }

	}

	report all

	set total [extract count.particles]
	set I     [expr [extract count.positions position=I]]
	set V     [expr [extract count.particles particle=V]]
	test tag=total.$c float=$total value=200 error=0
	test tag=I.$c     float=$I     value=100 error=0
	test tag=V.$c     float=$V     value=100 error=0

}
