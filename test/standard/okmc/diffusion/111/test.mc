param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"
param set type=int key=MC/General/snapshot.time.decade value=1

set size 40
set time 1e5

proc material { x y z } { return "S_Iron" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true
param set type=bool   key=MC/Mesh/periodic.y value=true
param set type=bool   key=MC/Mesh/periodic.z value=true


init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I10 coord={30 25 20}
insert defect=<111> ID=I100 coord={10 15 20}
insert defect=<111> ID=I1000 coord={5 10 15}

save lammps=nodist-defects

anneal time=$time temp=[expr 600 - 273.15] events=3e5

save lammps=nodist-defects append

test tag=time value=[expr 28.6*4] float=[extract time] error=0.01
