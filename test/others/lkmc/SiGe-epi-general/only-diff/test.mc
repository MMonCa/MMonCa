set T 600


set value(0) 6.10988
set value(1) 6.19598 
set value(2) 6.56812 
set value(3) 7.6034 
set value(4) 6.86634 
set value(5) 6.10988 
set value(6) 6.10988 
set value(7) 6.10988 


proc material { x y z } {
	if { $x <7 && $y > 3 && $y < 5 } { return "Gas" }
        if { $x >6 } { return "Silicon" }
        return "Gas"
} 


set sizeX [expr 24]
set sizeY [expr sqrt(2.)*0.5431*10]
set sizeZ [expr sqrt(2.)*0.5431*10]

param set type=map<string,string>   key=MC/General/materials value="Silicon Si Gas Gas"
param set type=map<string,bool>     key=Silicon/Models/defined value={ }
param set type=array<string,string> key=Silicon/Models/interactions value={ }
param set type=array<string>        key=Silicon/Models/interaction.result value={ }
param set type=bool                 key=Silicon/Models/epitaxy value=true
param set type=bool                 key=Silicon/Epitaxy/model.simplified  value=false
param set type=proc                 key=Silicon/Epitaxy/prefactor.precursor value={ {
return 0 
} }

#param set type=float                key=MC/General/snapshot.events        value=100000
#param set type=int                  key=MC/General/snapshot.time.decade   value=1


#proc snapshot { } { save lammps=nodist-Si append no.print }
init minx=0 miny=0 minz=0 maxx=$sizeX maxy=$sizeY maxz=$sizeZ material=material

#save lammps=nodist-Si
anneal time=1e15 temp=$T events=6000000
save lammps=nodist-Si

#exec gunzip -c Si-ref.dump.gz > nodist-Si-ref.dump
#lowmsg "Test PASSED -[exec diff nodist-Si.dump nodist-Si-ref.dump]-"

set roughness [lindex [extract ac.stdev] 0]

foreach i "0 1 2 3 4 5 6 7" {
	set depth [lindex [extract no.print ac.mean min.y=$i max.y=[expr $i+1]] 0]
	test tag=depth.$i float=$depth value=$value($i)0 error=0.02
}

test tag=rough float=$roughness value=0.58 error=0.03
