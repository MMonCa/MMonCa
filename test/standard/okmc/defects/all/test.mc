param set type=map<string,string>   key=MC/General/materials value="Silicon Si" 

set size 40
proc material { x y z } { return "Silicon" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<311>    ID=I38  coord={20 20 20}
insert defect=Void     ID=V39  coord={18 18 18}
insert defect=DLoop    ID=I40  coord={16 16 16}
insert defect=CCluster ID=C6I6 coord={14 14 14}

save lammps=nodist-defects
extract defects={ } file=nodist-defects.txt

test tag=Void     float=[extract count.particles defect=Void]     value=39 error=0
test tag=DLoop    float=[extract count.particles defect=DLoop]    value=40 error=0
test tag=<311>    float=[extract count.particles defect=<311>]    value=38 error=0
test tag=CCluster float=[extract count.particles defect=CCluster] value=6 error=0

anneal time=1e10 events=2e7 temp=700

save lammps=nodist-defects append
extract defects={ } file=nodist-defects2.txt
