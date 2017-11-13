param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe"

set size 143.5
set time 300

proc material { x y z } { return "S_Iron" }

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=true

param set type=float  key=S_Iron/Iron/I(capture.radius)    value=0.55 new

param set type=array<string,string> key=S_Iron/Models/interactions value={ {
	I+V         	0,0.55
	I+I         	ICluster,1.,0.55
	V+V         	VCluster,1.,0.55
	
	I+Gas       	true
	V+Gas       	true
    
	ICluster:*+I		true,0.55
	ICluster:*+V		true,0.55
	ICluster+ICluster	true,0.55
	ICluster+VCluster	true,0.55
	ICluster+<111>		true,0.55
	
	VCluster:*+I		true,0.55
	VCluster:*+V		true,0.55
	VCluster+VCluster	true,0.55
	VCluster+<111>      	true,0.55
	
	<111>:*+I		true,0.55
	<111>:*+V		true,0.55
	<111>+<111>         	true,0.55
} }
param set type=arrhenius key=S_Iron/Vacancy/V(migration)   value="5e-5 0.67"
param set type=arrhenius key=S_Iron/Iron/I(migration)      value="3.2e-3 0.34"

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

cascade file=electron.cascade format=B:C*.287:D*.287:E*.287 periodic fluence=4.8562e9 do.not.react

set X 0
set a 77.2
set b 1.030927

set FILE [open "results.txt" w]
close $FILE

set oldTotal 0
set oldC 0

save lammps=evolution
while { $X < 65 } {
      incr X
      set c [expr $a*$b]
#      lowmsg "Running for $c K"      
      set FILE [open "results.txt" a]
      anneal time=$time temp=[expr $c - 273.15]
      save lammps=evolution append no.print
      set total [extract count.particles no.print]
      set deriv [expr ($total - $oldTotal)/($c - $oldC)]
      set I     [extract count.particles particle=I defect=MobileParticle no.print]
      set V     [extract count.particles particle=V defect=MobileParticle no.print]
      set I2    [extract count.particles defect=ICluster ID=I2 no.print]
      set V2    [extract count.particles defect=VCluster ID=V2 no.print]
      set I3    [extract count.particles defect=ICluster ID=I3 no.print]
      set V3    [extract count.particles defect=VCluster ID=V3 no.print]
      set I111  [extract count.particles defect=<111> no.print]
      set oI    [expr [extract count.particles defect=ICluster no.print] -$I2 -$I3]
      set oV    [expr [extract count.particles defect=VCluster no.print] -$V2 -$V3]
      puts $FILE "$c $total $deriv $I $V $I2 $V2 $I3 $V3 $oI $oV $I111"
      lowmsg     "$c $total $deriv $I $V $I2 $V2 $I3 $V3 $oI $oV $I111"
      close $FILE
      set a $c
      set oldC $c
      set oldTotal $total
}

report reactions.interface
