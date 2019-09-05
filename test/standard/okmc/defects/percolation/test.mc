param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

set size 20
set time 300

proc material { x y z } {
	return "S_Iron"
}

#parameters
param set type=bool key=MC/Mesh/periodic.x value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster:*+He      value=true
param set type=array<string,string> key=S_Iron/Models/interactions index=He+He               value=HeCluster,1
param set type=array<string,string> key=S_Iron/Models/interactions index=HeCluster+HeCluster value=true
param set type=arrhenius key=S_Iron/Helium/He(migration) value={ 1 .4 }

proc snapshot { } { }

init minx=-20 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

save lammps=nodist-evolution


insert coord={ 10 10 10 } defect=HeCluster ID=He10
insert coord={ 10 11 10 } defect=HeCluster ID=He10

report all

set rad1 [extract defect.radius defect=HeCluster]
test tag=radius.1 float=0.269386 value=$rad1 error=.0001
test tag=1        float=2 value=[extract count.defects defect=HeCluster] error=0

save lammps=nodist-evolution
cascade file=cascade fluence=5e13 flux=1e8 voluminic periodic
report all
save lammps=nodist-evolution append

set rad2 [extract defect.radius defect=HeCluster]
test tag=radius.2 float=0.89 value=$rad2 error=.05
test tag=2        float=1 value=[extract count.defects defect=HeCluster] error=0
