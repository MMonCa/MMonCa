param set type=map<string,string>   key=MC/General/materials value="Iron Fe"
param set type=map<string,string>   key=Iron/Mechanics/eigenstrains value="V 0.5,0.00655379"
param set type=int   key=Mechanics/General/dimension value=3

param set type=array<string,string> key=Iron/Models/interactions index=<111>+Gas value=false
param set type=array<string,string> key=Iron/Models/interactions index=VCluster+Gas value=false

param set type=float key=MC/Mesh/spacing.x value=1.268
param set type=float key=MC/Mesh/spacing.y value=1.268
param set type=float key=MC/Mesh/spacing.z value=1.268

set elements 25

set T 550

proc material { x y z } { return "Iron" }
param set type=string key=Mechanics/General/model  value=Feliks
param set type=float  key=Mechanics/General/update.time.min value=1e-7
proc snapshot {} { }

init minx=0 miny=0 minz=0 maxx=[expr $elements*1.268] maxy=[expr $elements*1.268] maxz=[expr $elements*1.268] material=material

anneal time=1 temp=$T events=5


insert particle=V coord={ [expr $elements*1.268/2.] [expr $elements*1.268/2.] [expr $elements*1.268/2.] }
insert particle=V coord={ [expr $elements*1.268/3.] [expr $elements*1.268/3.] [expr $elements*1.268/3.] }
insert particle=V coord={ [expr $elements*1.268/1.5] [expr $elements*1.268/1.1] [expr $elements*1.268/1.1] }
insert particle=V coord={ [expr $elements*1.268/1.1] [expr $elements*1.268/1.5] [expr $elements*1.268/2.] }

anneal time=1 temp=$T events=5

extract stress name=xx file=nodist-sxx2.dat dimension=2

set FILE1 [exec cat sxx2.dat]
set FILE2 [exec cat nodist-sxx2.dat]

test arrays.2D one={ $FILE1 } two={ $FILE2 } error=.0009
