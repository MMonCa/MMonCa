param set type=map<string,string>   key=MC/General/materials value="Iron Fe Gas Gas"

proc material { x y z } {
	return "Iron"
}

param set type=bool key=MC/Mesh/periodic.x value=true
param set type=float key=MC/Mesh/spacing.x value=1
param set type=float key=MC/Mesh/spacing.y value=1
param set type=float key=MC/Mesh/spacing.z value=1
set alloyDensity [param get type=float key=Iron/Models/alloy.density]
set conc 0.5

proc alloy { x y z } {
    global alloyDensity
    global conc

    return [expr $conc * $alloyDensity]
}

param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I

init minx=0 miny=0 minz=0 maxx=10 maxy=10 maxz=10 material=material

profile name=Cr proc=alloy
insert particle=V coord={ 5 5 5 }

anneal time=1e25 temp=500 events=2e7

extract profile name=Cr dimension=3 file="gen.dat"
restart save=dump
