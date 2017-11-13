param set type=map<string,string>   key=MC/General/materials value="SiO2 SiO2 Nitride Ni"
set time 1

proc material { x y z } { 
	if { $x < 12 } { return "SiO2" }
	return "Nitride"
}

proc myName { x y z } { return 2e20 }

set size 30

param set type=arrhenius key=SiO2/Carbon/C(migration)    value="2.26e-3 0.791"
param set type=arrhenius key=SiO2/Carbon/C(formation)    value="4e25 0.0"
param set type=arrhenius key=Nitride/Carbon/C(formation) value="8e25 0.0"

init minx=0 miny=0 minz=0 maxx=24 maxy=$size maxz=$size material=material

profile name=C proc=myName

lowmsg [extract profile name=C]

anneal time=1e10 temp=800 events=2e8

lowmsg [extract profile.mobile name=C]
lowmsg [extract profile        name=C]

test float=[expr 2 * [extract count.particles particle=C material=SiO2]] value=[extract count.particles particle=C material=Nitride] error=0.10

test array="[extract profile name=C]"        value=1.3e20 error=0.40 init=0  end=11 tag=SiO2.sta
test array="[extract profile.mobile name=C]" value=1.3e20 error=0.10 init=0  end=11 tag=SiO2.dyn
test array="[extract profile name=C]"        value=2.6e20 error=0.25 init=13 end=24 tag=Nitride.sta
test array="[extract profile.mobile name=C]" value=2.6e20 error=0.10 init=13 end=24 tag=Nitride.dyn
