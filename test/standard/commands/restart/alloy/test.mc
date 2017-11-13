#param set type=map<string,string>   key=MC/General/materials value="Iron Fe"

#param set type=bool key=MC/Mesh/periodic.x value=true
#param set type=float key=MC/Mesh/spacing.x value=1
#param set type=float key=MC/Mesh/spacing.y value=1
#param set type=float key=MC/Mesh/spacing.z value=1

restart load=dump

extract profile name=Cr dimension=3 file=test.dat

set bPass yes
if { [catch { exec diff gen.dat test.dat} msg] } {
    set bPass no
}

test one=yes equal two=$bPass
