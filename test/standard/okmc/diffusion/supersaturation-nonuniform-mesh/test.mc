param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas" 
#Mon 25 Feb 2013 **** Test for supersaturation

set T 1000
set Bcon 5e20
set ev   3e7

set kB  8.6174e-5

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "Silicon"
}

proc spike { x y z } {       
       global Bcon
       return $Bcon
}

 

param set type=map<string,bool> key=Silicon/Models/particles index=I_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=I       value=true
param set type=map<string,int>  key=Silicon/Silicon/I(state.charge)        value="I 0"
param set type=map<string,bool> key=Silicon/Models/particles index=V_0,-,+ value=false
param set type=map<string,bool> key=Silicon/Models/particles index=V       value=true
param set type=map<string,int>  key=Silicon/Vacancy/V(state.charge)        value="V 0"
param set type=map<string,bool> key=Silicon/Models/particles index=BI_0,-    value=false
param set type=map<string,bool> key=Silicon/Models/particles index=BI        value=true
param set type=map<string,int>  key=Silicon/Boron/BI(state.charge)           value="BI -1"

param set type=array<string,string> key=Silicon/Models/interactions index=B+BI  value=false
param set type=array<string,string> key=Silicon/Models/interactions index=I+I value=false

param set type=arrhenius key=Silicon/Vacancy/V(formation) value={ 9.21e29 3.74 }
param set type=arrhenius key=Silicon/Silicon/I(formation) value={ 7.95e28 4 }
param set type=arrhenius key=Silicon/Vacancy/V(migration) value={ 5.00e-8 0.4 } new
param set type=arrhenius key=Silicon/Silicon/I(migration) value={ 2e-3 0.7 }    new
param set type=arrhenius key=Silicon/Boron/BI(formation)    value={ 3285 2.9 }
param set type=arrhenius key=Silicon/Boron/BI(migration)    value={ 1.7e-3 0.77 } new

set PfI   [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 0]
set EfI   [lindex [param get type=arrhenius key=Silicon/Silicon/I(formation)] 1]
set PfB   [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 0]
set EfB   [lindex [param get type=arrhenius key=Silicon/Boron/B(formation)] 1]

set Icon [expr $PfI * exp (- $EfI / ($kB * ($T + 273.15)))]

init linesx={-2 -1.5 -1 -0.5 0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} linesy={0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} linesz={0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} material=material

profile proc=spike name=B

anneal events=$ev temp=$T time=1

set B_diff [extract diffusivity macroscopic name=B material=Silicon]
lowmsg "Initial diffusion $B_diff"

# error was 0.05 in supersaturation
test tag=Icon.noSat array="[extract profile.mobile name=I material=Silicon]"  value=$Icon init=2 end=18 error=0.09

param set type=map<string,float> key=Gas_Silicon/Silicon/supersaturation.right index=I value=10

init linesx={-2 -1.5 -1 -0.5 0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} linesy={0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} linesz={0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 4.8 5.1 5.4 5.7 6} material=material

profile proc=spike name=B

anneal events=$ev temp=$T time=1

set I_10 [expr $Icon * 10]
set B_diff10 [expr $B_diff * 10]

# error was 0.05 in supersaturation
test tag=Icon.sat array="[extract profile.mobile name=I]" value=$I_10 init=2 end=18  error=0.08
test tag=diffB.sat float="[extract diffusivity macroscopic name=B material=Silicon]"  value=$B_diff10 error=0.30


