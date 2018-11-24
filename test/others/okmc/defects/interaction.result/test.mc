set size 20
param set type=map<string,string>   key=MC/General/materials value="Gas Gas Iron Fe"                

proc material { x y z } { return "Iron" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

param set type=string           key=Iron/<111>/from value=<111>
param set type=string           key=Iron/<111>/to   value=<111>
param set type=string           key=Iron/<100>/from value=<100>
param set type=string           key=Iron/<100>/to   value=<100>

param set type=array<string,string> key=Iron/Models/interactions value={ <111>+<111> true <100>+<100> true <111>+<100> true }
param set type=map<string,bool>     key=Iron/Models/defined value={ <100> true <111> true }
param set type=array<string>        key=Iron/Models/interaction.result value={ {
<100> + <100> = <100>
<111> + <111> = <111>

<100> ~=,.05 <111> = <111>
<100> <  <111> = <111>
<100> >  <111> = <100>

} }


init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<100> ID=I20 coord={ 15 15 15 }
insert defect=<100> ID=I25 coord={ 15 15 15 }

report all

test tag=100+100 float=[extract count.particles]            value=45 error=0
test tag=100+100 float=[extract count.defects defect=<100>] value=1  error=0
test tag=100+100 float=[extract count.defects defect=<111>] value=0  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I20 coord={ 15 15 15 }
insert defect=<111> ID=I25 coord={ 15 15 15 }

report all

test tag=111+111 float=[extract count.particles]            value=45 error=0
test tag=111+111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=111+111 float=[extract count.defects defect=<111>] value=1  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<100> ID=I20 coord={ 15 15 15 }
insert defect=<111> ID=I21 coord={ 15 15 15 }

report all

test tag=100~=111 float=[extract count.particles]            value=41 error=0
test tag=100~=111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=100~=111 float=[extract count.defects defect=<111>] value=1  error=0


#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<100> ID=I2 coord={ 15 15 15 }
insert defect=<111> ID=I21 coord={ 15 15 15 }

report all

test tag=100<111 float=[extract count.particles]            value=23 error=0
test tag=100<111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=100<111 float=[extract count.defects defect=<111>] value=1  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<100> ID=I21 coord={ 15 15 15 }
insert defect=<111> ID=I2  coord={ 15 15 15 }

report all

test tag=100>111 float=[extract count.particles]            value=23 error=0
test tag=100>111 float=[extract count.defects defect=<100>] value=1  error=0
test tag=100>111 float=[extract count.defects defect=<111>] value=0  error=0

#revert order just in case

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<100> ID=I25 coord={ 15 15 15 }
insert defect=<100> ID=I20 coord={ 15 15 15 }

report all

test tag=b100+100 float=[extract count.particles]            value=45 error=0
test tag=b100+100 float=[extract count.defects defect=<100>] value=1  error=0
test tag=b100+100 float=[extract count.defects defect=<111>] value=0  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I25 coord={ 15 15 15 }
insert defect=<111> ID=I20 coord={ 15 15 15 }

report all

test tag=b111+111 float=[extract count.particles]            value=45 error=0
test tag=b111+111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=b111+111 float=[extract count.defects defect=<111>] value=1  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I21 coord={ 15 15 15 }
insert defect=<100> ID=I20 coord={ 15 15 15 }

report all

test tag=b100~=111 float=[extract count.particles]            value=41 error=0
test tag=b100~=111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=b100~=111 float=[extract count.defects defect=<111>] value=1  error=0


#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I21 coord={ 15 15 15 }
insert defect=<100> ID=I2 coord={ 15 15 15 }

report all

test tag=b100<111 float=[extract count.particles]            value=23 error=0
test tag=b100<111 float=[extract count.defects defect=<100>] value=0  error=0
test tag=b100<111 float=[extract count.defects defect=<111>] value=1  error=0

#-------------
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

insert defect=<111> ID=I2  coord={ 15 15 15 }
insert defect=<100> ID=I21 coord={ 15 15 15 }

report all

test tag=b100>111 float=[extract count.particles]            value=23 error=0
test tag=b100>111 float=[extract count.defects defect=<100>] value=1  error=0
test tag=b100>111 float=[extract count.defects defect=<111>] value=0  error=0

