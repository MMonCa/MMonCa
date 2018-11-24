proc interpol { x x1 y1 x2 y2 } {
	set m [expr ($y2 - $y1) / ($x2 - $x1)]
	set c [expr ($y1 * $x2 - $y2 * $x1) / ($x2 - $x1)]
	return [expr $m * $x + $c]
}

proc material { x y z } {
	return "Iron"
}
param set type=map<string,string>   key=MC/General/materials value="Gas Gas Iron Fe"
param set type=bool                 key=Iron/Models/self.diffusion value=true

set conc 50
#SUPERLATTICE
proc alloy { x y z } {
       global conc
       global densityCm3
        if { $x < 8 } {
          return [expr $densityCm3 * ($conc + 0.05 * sin($x * 2 * 3.141593 / 16))]
        }
         return  [expr $densityCm3 * ($conc + 0.05 * sin($x * 2 * 3.141593 / 16))]
}

set sizeX  8
set sizeYZ 4
set T 500
set ev 4e6
set conc 0.5

param set type=bool   key=MC/Mesh/speed.up value=false

param set type=arrhenius key=Iron/Vacancy/V(formation) value={ 5.3e24 2.17 }
param set type=arrhenius key=Iron/Vacancy/V(migration) value={ 7e-2 0.7 }
param set type=arrhenius key=Iron/Iron/I(formation) value={ 0.01 10 }
param set type=arrhenius key=Iron/Iron/I(migration) value={ 0.01 10 }

param set type=float key=MC/Mesh/spacing.x value=1.0
param set type=float key=MC/Mesh/spacing.y value=1.0
param set type=float key=MC/Mesh/spacing.z value=1.0
param set type=bool key=MC/Mesh/periodic.x value=true
param set type=bool key=MC/Mesh/periodic.y value=true
param set type=bool key=MC/Mesh/periodic.z value=true

set densityCm3 [param get type=float key=Iron/Models/alloy.density]
set lambda [param get type=float key=Iron/Models/lambda]
set corrFact [expr 2. / $lambda]
set formP [lindex [param get type=arrhenius key=Iron/Vacancy/V(formation)] 0]
set formE [lindex [param get type=arrhenius key=Iron/Vacancy/V(formation)] 1]
set migP  [lindex [param get type=arrhenius key=Iron/Vacancy/V(migration)] 0]
set migE  [lindex [param get type=arrhenius key=Iron/Vacancy/V(migration)] 1]

param set type=float key=Iron/Vacancy/correlation.factor value=$corrFact
param set type=float key=Iron/Models/theta               value=2000

param set type=arrhenius key=Iron/Vacancy/alpha     value={ 5 0.1 }

proc alloy { x y z } {
       global conc
       global densityCm3
       return [expr $densityCm3 * ($conc + 0.05 * sin($x * 2 * 3.141593 / 8))]
}

init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

insert particle=V coord={ [expr $sizeX * 0.5] [expr $sizeYZ * 0.5] [expr $sizeYZ *0.5] }

profile name=Cr proc=alloy

extract mixing.enthalpy material=Iron file=Ho-nodist.mix

extract profile material=Iron file=Cr.init dimension=1 name=Cr
anneal time=1e25 temp=$T events=$ev
extract profile material=Iron file=Cr.fin dimension=1 name=Cr

test tag=alphaP array={ [extract profile name=Cr] } value=[expr 0.9 * $densityCm3] init=1 end=3 error=0.5
test tag=alpha array={ [extract profile name=Cr] } value=[expr 0.09 * $densityCm3] init=5 end=7 error=0.5


