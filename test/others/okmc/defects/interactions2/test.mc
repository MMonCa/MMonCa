param set type=map<string,string> key=MC/General/materials value="Iron Fe"

proc material { x y z } {
	return "Iron"
}

param set type=array<string,string> key=Iron/Models/interactions value={ {
	<111>+<111>	true
} }
param set type=map<string,bool> key=Iron/Models/defined value={ {
	ICluster true
	<111>   true
	<100>	true
} }

param set type=array<string> key=Iron/Models/interaction.result value={ {
	<111> == <111> = <111>
	<111> ~=,.05 <111> = <100>
	<111> > <111>  = <111>
	<111> < <111>  = <111>
} }

#----------------
init minx=0 miny=0 minz=0 maxx=15 maxy=15 maxz=15 material=material
insert defect=<111> ID=I200 coord={ 7 7 7 }
insert defect=<111> ID=I200 coord={ 7.1 7.1 7.1 }
anneal temp=-70 time=1e25 events=1
test tag="==" float=[extract count.defects defect=<111> ID=I400] value=1 error=0

#-----------------
init minx=0 miny=0 minz=0 maxx=15 maxy=15 maxz=15 material=material
insert defect=<111> ID=I200 coord={ 7 7 7 }
insert defect=<111> ID=I210 coord={ 7.1 7.1 7.1 }
anneal temp=-70 time=1e25 events=1
test tag="~=" float=[extract count.defects defect=<100> ID=I410] value=1 error=0

#-----------------
init minx=0 miny=0 minz=0 maxx=15 maxy=15 maxz=15 material=material
insert defect=<111> ID=I40 coord={ 7 7 7 }
insert defect=<111> ID=I20 coord={ 7.1 7.1 7.1 }
anneal temp=-70 time=1e25 events=1
test tag=">" float=[extract count.defects defect=<111> ID=I60] value=1 error=0

#-----------------
init minx=0 miny=0 minz=0 maxx=15 maxy=15 maxz=15 material=material
insert defect=<111> ID=I20 coord={ 7 7 7 }
insert defect=<111> ID=I40 coord={ 7.1 7.1 7.1 }
anneal temp=-70 time=1e25 events=1
test tag="<" float=[extract count.defects defect=<111> ID=I60] value=1 error=0
