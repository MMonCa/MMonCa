param set type=map<string,string>   key=MC/General/materials value="Silicon Si AmorphousSilicon aSi Gas Gas"

proc material { x y z } { 
	if { $x < 13.5 } { return "Gas" }
	return "Silicon"
}

set size 51

param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+B
param set type=array<string,string> key=Silicon/Models/interactions value=false index=B+BI
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+BI
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+V
param set type=array<string,string> key=Silicon/Models/interactions value=false index=V+BI
param set type=array<string,string> key=Silicon/Models/interactions value=false index=BI+I
param set type=array<string,string> key=Silicon/Models/interactions value=false index=I+BI

proc myName { x y z } { return 5e22 }

param set type=float	key=Silicon/Models/amorphization.threshold	value=6e22 new

init minx=0 miny=0 minz=0 maxx=24 maxy=$size maxz=$size material=material

lowmsg "Introducing damage below the threshold... This shouldn't amorphize."
profile name=I		proc=myName

report all

test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=27  tag=noAmorphI 


lowmsg "Allright, let's amorphize it!'"

profile name=V		proc=myName

report all

test array="[extract amorphous.fraction ]" value=1 error=0 init=13.5 end=27  tag=AmorphIV 
set thrs [param get type=float key=Silicon/Models/amorphization.threshold]
test array="[extract profile.damage ]" value=0 error=0 init=0 end=13.5 tag=profGas
lowmsg "threshold is $thrs"
test array="[extract profile.damage ]" value=$thrs error=0 init=13.5 end=27 tag=profSi

param set type=float	key=Silicon/Models/amorphization.threshold	value=1e22

init minx=0 miny=0 minz=0 maxx=27 maxy=$size maxz=$size material=material

lowmsg "What about changing the threshold?"

profile name=I		proc=myName

report all

test array="[extract amorphous.fraction]" value=0 error=0  init=0 end=13.5 tag=AmorphGasI 
test array="[extract amorphous.fraction]" value=1 error=0  init=13.5 end=27 tag=AmorphSiI 

init minx=0 miny=0 minz=0 maxx=27 maxy=$size maxz=$size material=material

lowmsg "These impurities should not amorphize, since no Is nor Vs are generated...."

profile name=B		proc=myName
profile name=BI		proc=myName

report all

test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=27  tag=noAmorphC

lowmsg "The last thing to check is to deactivate the amorphization model... "

param set type=float	key=Silicon/Models/amorphization.threshold	value=-1

init minx=0 miny=0 minz=0 maxx=27 maxy=$size maxz=$size material=material

profile name=I		proc=myName
profile name=V		proc=myName

report all

test array="[extract amorphous.fraction]" value=0 error=0 init=0 end=27  tag=noAmorphModel
test array="[extract profile.damage ]" value=1e23 error=0.05 init=13.5 end=27 tag=profSi



