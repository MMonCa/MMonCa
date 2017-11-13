lowmsg "This test checks if after SPER, the material is well updated. If it crashes, it may be due to bad updates in _state(LKMCAtom). Actually, it was created for that."
param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas" 

set maxx 12
set maxy 7.5
set maxz 7.5

set thres 2.97e20
proc material { x y z } {

	if {$x < 1.5} {
		return "Gas"
	}
	
		set res "Silicon"
		return $res
}

proc theProf { x y z } { 
	if {$x <= 7.5} {
		return [expr 18*3e20]
	}

}


set maxz    [expr sqrt(2.)*0.5431*10]
set maxy    [expr sqrt(2.)*0.5431*10]
param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true
set i 0 
set angle(0)  0
set radians "$angle($i).0*2.0*3.1415926535897931/360.0"
set S [expr sin($radians)]
set C [expr cos($radians)]
set R [expr sqrt(2.0)]
set waferorient "i $C  j [expr $S/$R] k [expr $S/$R]"
set flatorient  "i -$S j [expr $C/$R] k [expr $C/$R]"

param set type=map<string,float> key=Silicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=Silicon/Lattice/flat.orientation  value="$flatorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/wafer.orientation value="$waferorient"
param set type=map<string,float> key=AmorphousSilicon/Lattice/flat.orientation  value="$flatorient"


param set type=float		key=Silicon/Models/amorphization.threshold   value=$thres
init minx=0 miny=0 minz=0 maxx=$maxx maxy=$maxy maxz=$maxz material=material temp=-173

	profile name=I proc=theProf
lowmsg "This test checks if after SPER, the material is well updated. If it crashes, it may be due to bad updates in _state(LKMCAtom). Actually, it was created for that."
	test array="[extract amorphous.fraction]" value=1 error=0 init=2.25 end=6.75 tag=fractionBeforeSPER
	anneal time=500 temp=600
	extract amorphous.fraction file=profiles/after
lowmsg "This test checks if after SPER, the material is well updated. If it crashes, it may be due to bad updates in _state(LKMCAtom). Actually, it was created for that."
	test array="[extract amorphous.fraction]" value=0 error=0 init=0.75 end=11.25 tag=fractionAfterSPER



