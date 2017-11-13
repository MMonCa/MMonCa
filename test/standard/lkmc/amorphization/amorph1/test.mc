param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi" 

set maxx 9
set maxy 9
set maxz 9

set type 1
proc material { x y z } {

		set res "Silicon"

}

proc theProf { x y z } { 
	return [expr 1e25]
}
param set type=bool		key=MC/Mesh/periodic.x value=true
param set type=float		key=MC/Mesh/spacing.x  value=1.5
param set type=float		key=MC/Mesh/spacing.y  value=1.5
param set type=float		key=MC/Mesh/spacing.z  value=1.5
param set type=float		key=Silicon/Models/amorphization.threshold   value=$type

init minx=0 miny=0 minz=0 maxx=$maxx maxy=$maxy maxz=$maxz material=material temp=-173
	lowmsg "testing some atoms post-amorphization"
	save lammps="nodist-Si" scale=10 no.print
	insert particle=I coord={2.5 2.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print


	test float="[lindex [extract ac.max] 0]" value=4.61635 error=0 tag=maxX1
	test float="[lindex [extract ac.max] 1]" value=4.60836 error=0 tag=maxY1
	test float="[lindex [extract ac.max] 2]" value=4.60836 error=0 tag=maxZ1

	test float="[lindex [extract ac.min] 0]" value=1.49353 error=0 tag=minX1
	test float="[lindex [extract ac.min] 1]" value=1.3441 error=0 tag=minY1
	test float="[lindex [extract ac.min] 2]" value=1.3441 error=0 tag=minZ1

	test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=1.5 tag=noAmorphBegin1
	test array="[extract amorphous.fraction]" value=0.111111 error=0 init=2.25 end=3.75 tag=Amorph1
	test array="[extract amorphous.fraction]" value=0 error=0 init=4.5 end=9 tag=noAmorphEnd1

	lowmsg "testing some atoms removal"

	insert particle=I coord={2.5 2.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={2.5 2.5 6}
	save lammps="nodist-Si" scale=10 append no.print

	test float="[lindex [extract ac.max] 0]" value=4.61635 error=0 tag=maxX2
	test float="[lindex [extract ac.max] 1]" value=4.60836 error=0 tag=maxY2
	test float="[lindex [extract ac.max] 2]" value=8.83268 error=0 tag=maxZ2

	test float="[lindex [extract ac.min] 0]" value=1.49353 error=0 tag=minX2
	test float="[lindex [extract ac.min] 1]" value=1.3441 error=0 tag=minY2
	test float="[lindex [extract ac.min] 2]" value=0 error=0 tag=minZ2

	test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=1.5 tag=noAmorphBegin2
	test array="[extract amorphous.fraction]" value=0.333333 error=0 init=2.25 end=3.75 tag=Amorph2
	test array="[extract amorphous.fraction]" value=0 error=0 init=4.5 end=9 tag=noAmorphEnd2

	insert particle=I coord={2.5 0.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={2.5 0.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={2.5 0.5 6}
	save lammps="nodist-Si" scale=10 append no.print

	insert particle=I coord={2.5 6 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={2.5 6 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={2.5 6 6}
	save lammps="nodist-Si" scale=10 append no.print

	set hola [extract ac.max]
	set adios [extract ac.min]
	lowmsg "hola $hola"
	lowmsg "adios $adios"

	lowmsg "checking interface"
	test float="[lindex [extract ac.max max.x=3] 0]" value=1.49353 error=0 tag=maxX3

	test float="[lindex [extract ac.min min.x=2] 0]" value=4.61635 error=0 tag=minX3

	test float="[lindex [extract ac.max max.x=4.4 min.x=1.76155] 0]" value=1.76155 error=0 tag=noAtoms
#this absurd value confirms that there are no atoms between max.x and min.x


	test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=1.5 tag=noAmorphBegin3
	test array="[extract amorphous.fraction]" value=1 error=0 init=2.25 end=3.75 tag=Amorph3
	test array="[extract amorphous.fraction]" value=0 error=0 init=4.5 end=9 tag=noAmorphEnd3


	lowmsg "testing some more"
	insert particle=I coord={0.5 2.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 2.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 2.5 6}
	save lammps="nodist-Si" scale=10 append no.print


	insert particle=I coord={0.5 0.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 0.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 0.5 6}
	save lammps="nodist-Si" scale=10 append no.print

	insert particle=I coord={0.5 6 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 6 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={0.5 6 6}
	save lammps="nodist-Si" scale=10 append no.print

	insert particle=I coord={6 2.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={6 2.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={6 2.5 6}
	save lammps="nodist-Si" scale=10 append no.print

	insert particle=I coord={6 0.5 2.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={6 0.5 0.5}
	save lammps="nodist-Si" scale=10 append no.print
	insert particle=I coord={6 0.5 6}
	save lammps="nodist-Si" scale=10 append no.print

	test float="[lindex [extract ac.max] 0]" value=7.46763 error=0 tag=maxXFin
	test float="[lindex [extract ac.max] 1]" value=7.48858 error=0 tag=maxYFin
	test float="[lindex [extract ac.max] 2]" value=8.83268 error=0 tag=maxZFin

	test float="[lindex [extract ac.min] 0]" value=4.61635 error=0 tag=minXFin
	test float="[lindex [extract ac.min] 1]" value=4.60836 error=0 tag=minYFin
	test float="[lindex [extract ac.min] 2]" value=0 error=0 tag=minZFin

	test array="[extract amorphous.fraction]" value=1 error=0 init=0 end=4.4 tag=noAmorphBeginFin
	test array="[extract amorphous.fraction]" value=0.666667 error=0 init=4.5 end=7.5 tag=AmorphFin
	test array="[extract amorphous.fraction]" value=1 error=0 init=7.6 end=9 tag=noAmorphEndFin

init minx=0 miny=0 minz=0 maxx=$maxx maxy=$maxy maxz=$maxz material=material temp=-173
	lowmsg "this profile should not let any atom behind"
	profile name=I proc=theProf

	test array="[extract amorphous.fraction]" value=1 error=0 init=0 end=10 tag=amorphProfile
	test float="[lindex [extract ac.max] 0]" value=0 error=0 tag=maxXProf
	test float="[lindex [extract ac.max] 1]" value=0 error=0 tag=maxYProf
	test float="[lindex [extract ac.max] 2]" value=0 error=0 tag=maxZProf
	test float="[lindex [extract ac.min] 0]" value=9 error=0 tag=minXProf
	test float="[lindex [extract ac.min] 1]" value=9 error=0 tag=minYProf
	test float="[lindex [extract ac.min] 2]" value=9 error=0 tag=minZProf




