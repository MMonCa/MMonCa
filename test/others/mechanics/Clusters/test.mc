param set type=map<string,string>   key=MC/General/materials value="Iron Fe"
#param set type=map<string,string>   key=Iron/Mechanics/eigenstrains value={ {
#V      0.5,0.00655379
#He2V   0.5284190389705387,-0.0038966613396029626
#He2V2  0.4235275503426886,0.008050695320718804
#} }

param set type=map<string,string>   key=Iron/Mechanics/eigenstrains value={ {
V      0.5,0.0065
He2V   0.4999,0.0065
He2V2  0.5,-0.0060
} }

param set type=array<string,string> key=Iron/Models/interactions index=<111>+Gas value=false
param set type=array<string,string> key=Iron/Models/interactions index=VCluster+Gas value=false
param set type=map<string,bool>     key=Iron/Models/defined      index=HeCluster    value=true
param set type=map<string,bool>     key=Iron/Models/particles    index=He           value=true
param set type=map<string,bool>     key=Iron/Models/particles    index=HeV          value=true

param set type=int   key=Mechanics/General/dimension value=3

param set type=float key=MC/Mesh/spacing.x value=1.268
param set type=float key=MC/Mesh/spacing.y value=1.268
param set type=float key=MC/Mesh/spacing.z value=1.268

set elements 15

set T 550

proc material { x y z } { return "Iron" }
param set type=string key=Mechanics/General/model  value=Feliks
param set type=float  key=Mechanics/General/update.time.min value=1e-9

proc snapshot {} { }
set l [expr $elements*1.268]

init minx=0 miny=0 minz=0 maxx=[expr $l] maxy=[expr $l] maxz=[expr $l] material=material

anneal time=1 temp=$T events=5

insert               particle=V  coord={ [expr $l*.25] [expr $l*.25] [expr $l/2.] }
insert defect=HeCluster ID=He2V  coord={ [expr $l*.75] [expr $l*.25] [expr $l/2.] }
insert defect=HeCluster ID=He2V2 coord={ [expr $l*.25] [expr $l*.75] [expr $l/2.] }
insert               particle=V  coord={ [expr $l*.75] [expr $l*.75] [expr $l/2.] }

anneal time=1 temp=$T events=1

extract stress name=xx file=nodist-sxx2.dat dimension=2

set FILE1 [exec cat sxx2.dat]
set FILE2 [exec cat nodist-sxx2.dat]

test arrays.2D one={ $FILE1 } two={ $FILE2 } error=.0009
