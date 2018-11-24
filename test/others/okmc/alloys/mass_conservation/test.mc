param set type=map<string,string>  key=MC/General/materials  value="Iron Fe Gas Gas"

proc material { x y z } {
	if { $x < 0 } {
		return "Gas"
	}
	return "Iron"
}

proc alloy { x y z } {
	return 1e22
}

set size 18

param set type=bool   key=Iron/Models/self.diffusion value=true

init minx=-1.5 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material	
profile proc=alloy name=Cr

anneal time=1e25 temp=500 events=1e7


test tag=Fe.balance array={ [extract balance name=A.atoms material=Iron file=Fe.bal] } value=0 init=0 end=18 error=0
test tag=Cr.balance array={ [extract balance name=B.atoms material=Iron file=Cr.bal] } value=0 init=0 end=18 error=0


