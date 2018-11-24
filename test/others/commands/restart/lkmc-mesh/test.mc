set otime1  0.0100012
set odepth1 5.58011
set orough1 1.44814

set otime2  0.0200058
set odepth2 5.58336
set orough2 1.42926

restart load=dump
report all

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

test tag=time1  float=$time  value=$otime1  error=0
test tag=depth1 float=$depth value=$odepth1 error=0
test tag=rough1 float=$rough value=$orough1 error=0

anneal temp=600 time=0.01 epitaxy="Si 1."

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

save lammps=nodist-epi

test tag=time2  float=$time  value=$otime2  error=0.02
test tag=depth2 float=$depth value=$odepth2 error=0.07
test tag=rough2 float=$rough value=$orough2 error=0.07
