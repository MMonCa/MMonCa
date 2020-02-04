set results(0.05,1600) 1263.753538935439
set results(0.05,1400) 54.6540940507993
set results(0.05,800) 0.00017910551980573506
set results(0.24,1600) 1991.9238100144241
set results(0.24,1400) 205.57011656188945
set results(0.24,800) 0.0006813938120349445

param set type=map<string,string>   key=MC/General/materials value="GalliumArsenide GaAs Gas Gas"
param set type=map<string,string>   key=MC/Particles/elements value="As Arsenic,33,74.9216 Ga Gallium,11,22"
param set type=array<string,string> key=GalliumArsenide/Models/interactions value={ }
param set type=map<string,bool>     key=GalliumArsenide/Models/particles value={As false}
param set type=bool                 key=GalliumArsenide/Models/epitaxy value=true

set FILE [open "nodist-RP-CVD5.dat" w]
close $FILE

set sizeZ    		[expr sqrt(2.)*.5431*26]
set sizeY    		[expr sqrt(2.)*.5431*26]
set interface		52
set growth		50

proc material { x y z } {
	global interface
	if { $x < $interface } { 
		set res "Gas" 
	} else {
		set res "GalliumArsenide"
	}
	return $res
}

set Vepi	1.5E10
set Eepi	4.0
set Vdes	1.0E6
set Edes	0.32
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/prefactor.epi	 index=Ga value=$Vepi
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/barrier.epi		 index=Ga value=$Eepi
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/prefactor.desorption index=Ga value=$Vdes
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/barrier.desorption   index=Ga value=$Edes

param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/prefactor.epi	 index=As value=$Vepi
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/barrier.epi		 index=As value=$Eepi
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/prefactor.desorption index=As value=$Vdes
param set type=map<string,float> 	key=GalliumArsenide/Epitaxy/barrier.desorption   index=As value=$Edes

set tempList	[list 1600 1400 800]
set presList	[list 0.05 0.24]

foreach P $presList {
foreach T $tempList {

	init minx=-2 miny=0 minz=0 maxx=54 maxy=$sizeY maxz=$sizeZ material=material
	save lammps=nodist-$P-$T

	set grown [expr $interface - 0.5]
	anneal time=0.01 temp=$T depth=$grown epitaxy="Ga $P As $P"
	save lammps=nodist-$P-$T append

	set ini_time [extract time]
	set ini_depth(100) [lindex [extract ac.mean] 0]


	anneal time=0.01 temp=$T depth=$growth epitaxy="Ga $P As $P"
	set end_time [extract time]
	set elapsed_time [expr $end_time - $ini_time]
	set end_depth(100) [lindex [extract ac.mean] 0]
	set GR(100) [expr abs($ini_depth(100)-$end_depth(100))/$elapsed_time*60]

	set FILE [open "nodist-RP-CVD5.dat" a]
	puts $FILE "set results($P,$T) $GR(100)"
	close $FILE

	save lammps=nodist-$P-$T append
	test tag=$P,$T float=$GR(100) value=$results($P,$T) error=0.1
} }

