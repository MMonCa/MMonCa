param set type=map<string,string>   key=MC/General/materials value="S_Iron Fe Gas Gas"
param set type=int key=MC/General/domains value=2

set size 50
set time 300

proc material { x y z } {
	if { $x < 0 } { return "Gas" }
	return "S_Iron"
}

#parameters
param set type=bool   key=MC/Mesh/periodic.x value=false

proc snapshots { } {
 report events
}

param set type=proc   key=S_Iron/<111>/migration             value={ {
        set list "" 
        for { set size 2 } { $size < 500 } { incr size } {
                set pref [expr 3.5e-4 + 1.7e-3/pow($size,1.7)]
                set ener 1.0
                if { $size == 2 } { set pref 8.2e-3; set ener 0.42 }
                if { $size == 3 } { set pref 8.2e-3; set ener 0.43 }
                if { $size == 4 } { set pref 8.2e-3; set ener 0.43 }
                lappend list I$size
                lappend list $pref
                lappend list $ener
        }
} }


init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

cascade file=data/cascades periodic fluence=5e14 flux=1e12 temp=-150

set X 0
set a 77.2
set b 1.030927

set FILE [open "nodist-results.txt" w]
close $FILE

set oldTotal 0
set oldC 0

while { $X < 30 } {
      incr X
      set c [expr $a*$b]
      lowmsg "Running for $c K"
      set FILE [open "nodist-results.txt" a]
      cascade file=data/cascades periodic fluence=1e13 flux=1e12 temp=-150

      anneal time=$time temp=[expr $c - 273.15]
      set total [extract count.particles]
      set deriv [expr ($total - $oldTotal)/($c - $oldC)]
      set I     [extract count.particles particle=I defect=MobileParticle]
      set V     [extract count.particles particle=V defect=MobileParticle]
      set I2    [extract count.particles defect=ICluster ID=I2]
      set V2    [extract count.particles defect=VCluster ID=V2]
      set I3    [extract count.particles defect=ICluster ID=I3]
      set V3    [extract count.particles defect=VCluster ID=V3]
      set oI    [expr [extract count.particles defect=ICluster] -$I2 -$I3]
      set oV    [expr [extract count.particles defect=VCluster] -$V2 -$V3]
      puts $FILE "$c $total $deriv $I $V $I2 $V2 $I3 $V3 $oI $oV 0"
      lowmsg     "$c $total $deriv $I $V $I2 $V2 $I3 $V3 $oI $oV 0V"
      close $FILE
      set a $c
      set oldC $c
      set oldTotal $total
}
set Is [extract count.particles particle=I]
set Vs [extract count.particles particle=V defect=MobileParticle]

test tag=count.I float=$Is value=0    error=0.1
test tag=count.V float=$Vs value=7803 error=0.12
