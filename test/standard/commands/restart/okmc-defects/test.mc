restart load=dump

report defects

test array="[extract profile name=I]" value=4e20 error=0.06  init=1 end=11 tag=I.par
test array="[extract profile name=C]" value=8e20 error=0.06  init=1 end=11 tag=C.par
test array="[extract profile name=I]" value=0 error=0.  init=13 end=24 tag=I.par
test array="[extract profile name=C]" value=0 error=0.  init=13 end=24 tag=C.par

set size 50
test float=[extract count.defects ID=C2 defect=CCluster] value=[expr 12*$size*$size*1e-21*4e20] error=0.05 tag=Defects
test float=[extract count.particles particle=B] value=1 error=0 tag=B.interface
