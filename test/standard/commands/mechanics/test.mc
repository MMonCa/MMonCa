param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas"
param set type=float key=AmorphousSilicon/Mechanics/amorphous.expansion value=0.02

set T 550

proc material { x y z } {
        set res "Unknown"
	if { $x < 0 } { 
		set res "Gas" 
	} elseif { $x > 30 && $x < 50 && $y > 50 && $y < 150 } { 
		set res "AmorphousSilicon" } else {
		set res "Silicon"
	}
	return $res
}

param set type=bool   key=MC/Mesh/periodic.y       value=false
param set type=bool   key=MC/Mesh/periodic.z       value=true
param set type=string key=Mechanics/General/model  value=Feliks

proc snapshot {} { }

init minx=31 miny=0 minz=0 maxx=80 maxy=200 maxz=30 material=material

anneal time=1 temp=$T events=5

extract stress name=xy file=nodist-sxy.dat dimension=2
extract stress name=xx file=nodist-sxx.dat dimension=2
extract stress name=yy file=nodist-syy.dat dimension=2
extract strain name=xy file=nodist-xy.dat dimension=2
extract strain name=xx file=nodist-xx.dat dimension=2
extract strain name=yy file=nodist-yy.dat dimension=2

set FILE1 [exec cat data/sxy.ref]
set FILE2 [exec cat nodist-sxy.dat]

test arrays.2D one={ $FILE1 } two={ $FILE2 } error=.003
