param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=bool                 key=Silicon/Models/self.diffusion  value=true
param set type=int                key=MC/General/random.seed value=[clock seconds]

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}

set alloyConcCm3 2e22
proc alloy { x y z } {
	global alloyConcCm3
    if { $x > 15 && $x < 17 } { return $alloyConcCm3; }
    return 0;
}

set sizeX  34.5
set sizeYZ 9
set T 1000
set time 200
set Tot [param get type=float key=Silicon/Models/atomic.density]
set kB 8.62e-5
set toK 273.15
set xo 16.0
set comp [expr $alloyConcCm3 / $Tot ]

param set type=float     key=MC/Mesh/spacing.x value=1

param set type=map<string,bool> key=Silicon/Models/particles index=I_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=I       value=true
param set type=map<string,int>  key=Silicon/Silicon/I(state.charge)   value="I 0"
param set type=map<string,bool> key=Silicon/Models/particles index=V_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=V       value=true
param set type=map<string,int>  key=Silicon/Vacancy/V(state.charge)        value="V 0"
param set type=array<string,string> key=Silicon/Models/interactions index=C+I value=false
param set type=array<string,string> key=Silicon/Models/interactions index=C+V value=false

param set type=arrhenius key=Silicon/Silicon/I(migration) value={ 2e-3 0.7 } new
param set type=arrhenius key=Silicon/Silicon/I(formation) value={ 7.95e27 4 } new
param set type=arrhenius key=Silicon/Silicon/alpha value={ 1 0 }
param set type=arrhenius key=Silicon/Vacancy/V(migration) value={ 5.e-8 0.4 } new
param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 9.2e29 3.74 } new
param set type=arrhenius key=Silicon/Vacancy/alpha value={ 1 0 }

set fI      [param get type=float key=Silicon/Silicon/correlation.factor]
set CI_pref [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 0]
set CI_ener [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set DI_pref [lindex [param get type=arrhenius key=Silicon/Silicon/I(migration)] 0]
set DI_ener [lindex [param get type=arrhenius key=Silicon/Silicon/I(migration)] 1]
set aI_pref [lindex [param get type=arrhenius key=Silicon/Silicon/alpha] 0]
set aI_ener [lindex [param get type=arrhenius key=Silicon/Silicon/alpha] 1]

set fV      [param get type=float key=Silicon/Vacancy/correlation.factor]
set CV_pref [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 0]
set CV_ener [lindex [param get type=arrhenius key=Silicon/Vacancy/V(formation)] 1]
set DV_pref [lindex [param get type=arrhenius key=Silicon/Vacancy/V(migration)] 0]
set DV_ener [lindex [param get type=arrhenius key=Silicon/Vacancy/V(migration)] 1]
set aV_pref [lindex [param get type=arrhenius key=Silicon/Vacancy/alpha] 0]
set aV_ener [lindex [param get type=arrhenius key=Silicon/Vacancy/alpha] 1]

set CI [expr $CI_pref * exp(-$CI_ener/($kB*($T+$toK)))]
set DI [expr $DI_pref * exp(-$DI_ener/($kB*($T+$toK)))]
set aI [expr $aI_pref * exp(-$aI_ener/($kB*($T+$toK)))]
set CV [expr $CV_pref * exp(-$CV_ener/($kB*($T+$toK)))]
set DV [expr $DV_pref * exp(-$DV_ener/($kB*($T+$toK)))]
set aV [expr $aV_pref * exp(-$aV_ener/($kB*($T+$toK)))]

set Diff_GeI [expr $aI * $fI * $DI * $CI / $Tot / (1 - $comp + $aI * $comp)]
set Diff_GeV [expr $aV * $fV * $DV * $CV / $Tot / (1 - $comp + $aV * $comp)]
set Diff_Ge  [expr $Diff_GeI + $Diff_GeV]

set Dt       [expr $Diff_Ge * $time * 1e14]

param set type=array<string,string> key=Silicon/Models/interactions value={ {
	I+Gas true
	V+Gas true
} }


init minx=-1.5 miny=0 minz=0 maxx=$sizeX maxy=$sizeYZ maxz=$sizeYZ material=material

profile name=Ge  proc=alloy

extract profile name=Ge file=Ge.init

anneal time=$time temp=$T

extract profile name=Ge file=Ge-nodist.profile

array set Ger [extract profile name=Ge]

set area 0
set pi 3.14159265
foreach {key value} [array get Ger]  {
    set area  [expr $area + $value]
}

set F [open "Ge-nodist.th" "w"]
for {set i 0} {$i < [expr $sizeX * 10]} {incr i} {
    set Gauss($i)  [expr $area/(2*sqrt($pi*$Dt)) * exp(0 - ($i/10.0 - $xo)**2/(4*$Dt)) ]
    puts $F "[expr $i / 10.] $Gauss($i)"
}
close $F
foreach {i} {130 140 150 160 170 180 190} {
	set j [expr int($i / 10.)]
	test tag=SelfDiffusion.in.$i float=$Ger($j) value=$Gauss($i) error=1
}
