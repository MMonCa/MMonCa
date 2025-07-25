param set type=map<string,string>   key=MC/General/materials value="Silicon Si SiO2 SiO2"

param set type=int key=MC/General/domains value=2

proc material { x y z } { 
  if { ($x - 3) * ($x - 3) + $y * $y + $z * $z < 9 } { return "SiO2"; } 
  if { $x < 1.5 } { return "SiO2" }
  return "Silicon"
}

init minx=0 miny=0 minz=0 maxx=6 maxy=6 maxz=6 material=material

extract material.map filename="uniform-2.vtk"
test one="SiO2" equal two=[extract material.location x=0.1 y=0.1 z=0.1] tag=uniform_SiO2
test one="SiO2" equal two=[extract material.location x=0.1 y=5.0 z=1.1] tag=uniform_SiO2
test one="Silicon" equal two=[extract material.location x=5.1 y=5.0 z=3.1] tag=uniform_Si



init mesh="mesh.json"

extract material.map filename="nonuniform-json-2.vtk"
test one="Silicon" equal two=[extract material.location x=0.1 y=0.1 z=0.1] tag=nonuniform_json_SiO2
test one="Silicon" equal two=[extract material.location x=0.1 y=5.0 z=1.1] tag=nonuniform_json_SiO2
test one="SiO2" equal two=[extract material.location x=5.1 y=5.0 z=3.1] tag=nonuniform_json_Si



init linesx={0 1 2 3 4 5 6} linesy={0 0.5 1 1.5 2 3 6} linesz={0 0.25 0.75 1.75 3.75 6}  material=material

extract material.map filename="nonuniform-2.vtk"
test one="SiO2" equal two=[extract material.location x=0.1 y=0.1 z=0.1] tag=nonuniform_SiO2
test one="SiO2" equal two=[extract material.location x=0.1 y=5.0 z=1.1] tag=nonuniform_SiO2
test one="Silicon" equal two=[extract material.location x=5.1 y=5.0 z=3.1] tag=nonuniform_Si
