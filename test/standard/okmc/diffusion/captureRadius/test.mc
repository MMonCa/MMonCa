param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe" 
proc material { x y z } { return "S_Iron" }
set size  10

param set type=array<string,string> key=S_Iron/Models/interactions value=true,.19 index=C+I
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
insert coord={ 5 5.0 5 } particle=C
insert coord={ 5 5.2 5 } particle=I
test tag=no float={ [extract count.particles particle=CI] } value=0 error=0


param set type=array<string,string> key=S_Iron/Models/interactions value=true,.21 index=C+I
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material  
insert coord={ 5 5.0 5 } particle=C                                                                 
insert coord={ 5 5.2 5 } particle=I
test tag=yes float={ [extract count.particles particle=CI] } value=1 error=0                                   
