set time 1

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}

proc Geprof { x y z } {
	if { $x < 30 } { return 1e21; }
	return 0;
}

set alloyDen [param get type=float key=Silicon/Models/alloy.density]
proc Crprof { x y z } {
	global alloyDen;
	if { $x < 30 } {
		return $alloyDen;
	} else {
		return 0;
	}
}

set sizeX  60
set sizeYZ 12

init material=material minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ

profile name=Ge proc=Geprof

test tag="Ge.profile" array={ [extract profile name=Ge] } value=1e21 init=1  end=29 error=.06
test tag="Ge.profile" array={ [extract profile name=Ge] } value=0    init=31 end=59 error=.0

param set type=string    key=Silicon/Models/alloy value=Chromium

init material=material minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ

profile name=Cr proc=Crprof

set cellSide [param get type=float key=MC/Mesh/spacing.x]
set nAtm [expr $alloyDen * $cellSide ** 3 * 1e-21]

test tag="Cr.B.atoms" array={ [extract profile name=B.atoms] } value=$nAtm init=1 end=29 error=.01
test tag="Cr.A.atoms" array={ [extract profile name=A.atoms] } value=$nAtm init=31 end=59 error=.01



