param set type=map<string,string>   key=MC/General/materials value="S_Iron SFe"
proc material { x y z } { return "S_Iron" }
set size 20

init minx=0 miny=0 minz=0 maxx=$size maxy=$size maxz=$size material=material

