param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"

proc tst { tag str1 str2 } {
     if { $str1 == $str2 } { 
     	lowmsg "$tag $str1 -- $str2 OK"
     } else {
       lowmsg "$tag $str1 != $str2 FAILED"
       error "$tag is $str1 instead of $str2 TEST FAILED"
     }
}

proc material { x y z } { return "S_Iron" }
set size 20

test tag=first value=1 float=1 error=0

param set type=map<string,bool>     key=S_Iron/Models/defined     value=true index=HeCluster
param set type=array<string,string> key=S_Iron/Models/interactions value={ }

# clusters are He2 He2V HeV2 He2V2

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
tst a [param get.reaction material=S_Iron index=HeCluster:He2+V] false
tst b [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst c [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false
tst d [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false
tst f [param get.reaction material=S_Iron index=HeCluster+VCluster] false
tst g [param get.reaction material=S_Iron index=HeCluster+HeCluster] false

param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster:He*+V
param set type=array<string,string> key=S_Iron/Models/interactions value=false index=HeCluster:He*V*+V
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
tst h [param get.reaction material=S_Iron index=HeCluster:He2+V] true
tst i [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst j [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false
tst k [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false
tst l [param get.reaction material=S_Iron index=HeCluster+VCluster] false
tst m [param get.reaction material=S_Iron index=HeCluster+HeCluster] false

param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster:He*+V
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
tst n [param get.reaction material=S_Iron index=HeCluster:He2+V] true
tst o [param get.reaction material=S_Iron index=HeCluster:He2V+V] true
tst p [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false
tst q [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false
tst r [param get.reaction material=S_Iron index=HeCluster+VCluster] false
tst s [param get.reaction material=S_Iron index=HeCluster+HeCluster] false

param set type=array<string,string> key=S_Iron/Models/interactions value=false index=HeCluster:He*+V
param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster:HeV2+V
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material  
tst t [param get.reaction material=S_Iron index=HeCluster:He2+V] false
tst u [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst v [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false
tst w [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false 
tst x [param get.reaction material=S_Iron index=HeCluster+VCluster] false
tst y [param get.reaction material=S_Iron index=HeCluster+HeCluster] false

param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster+HeCluster
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material
tst 0 [param get.reaction material=S_Iron index=HeCluster:He2+V] false
tst 1 [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst 2 [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false 
tst 3 [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false 
tst 4 [param get.reaction material=S_Iron index=HeCluster+VCluster] false
tst 5 [param get.reaction material=S_Iron index=HeCluster+HeCluster] true

param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster+VCluster
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material  
tst 6 [param get.reaction material=S_Iron index=HeCluster:He2+V] false
tst 7 [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst 8 [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false 
tst 9 [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false
tst A [param get.reaction material=S_Iron index=HeCluster+VCluster] true
tst B [param get.reaction material=S_Iron index=HeCluster+HeCluster] true 
tst C [param get.reaction material=S_Iron index=HeCluster+Gas] false

param set type=array<string,string> key=S_Iron/Models/interactions value=true  index=HeCluster+Gas
init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material  
tst D [param get.reaction material=S_Iron index=HeCluster:He2+V] false
tst E [param get.reaction material=S_Iron index=HeCluster:He2V+V] false
tst F [param get.reaction material=S_Iron index=HeCluster:He2V2+V] false 
tst G [param get.reaction material=S_Iron index=HeCluster:HeV2+V] false
tst H [param get.reaction material=S_Iron index=HeCluster+VCluster] true
tst I [param get.reaction material=S_Iron index=HeCluster+HeCluster] true 
tst J [param get.reaction material=S_Iron index=HeCluster+Gas] true
