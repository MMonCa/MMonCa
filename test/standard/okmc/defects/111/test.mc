param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"


set size 40
set time 1e6

proc material { x y z } { return "S_Iron" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false
param set type=bool   key=MC/Mesh/periodic.y value=false
param set type=bool   key=MC/Mesh/periodic.z value=false


init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I20 coord={30 25 20}
insert defect=VCluster ID=V20 coord={18 18 18}

save lammps=nodist-defects

anneal time=$time temp=[expr 600 - 273.15] events=3e6

test tag=perpendicular float=[extract count.defects defect=VCluster] value=1 error=0

extract defects={ } file=nodist-out

#####################################################

param set type=string key=S_Iron/<111>/migration.type value=3d

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I20 coord={30 25 20}
insert defect=VCluster ID=V20 coord={18 18 18}

anneal time=$time temp=[expr 600 - 273.15] events=3e6

test tag=3D float=[extract count.defects defect=VCluster] value=0 error=0

