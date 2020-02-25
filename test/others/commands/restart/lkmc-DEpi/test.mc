
set otime1  1.00002
set odepth1 -8.20179
set orough1 0.0930396

set otime2  1.50002
set odepth2 -15.0963
set orough2 0.0608406

restart load=dump
report all

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

test tag=time1  float=$time  value=$otime1  error=0
test tag=depth1 float=$depth value=$odepth1 error=0
test tag=rough1 float=$rough value=$orough1 error=0

anneal temp=850 time=.5 epitaxy="Si 1."

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

save lammps=nodist-epi

test tag=time2  float=$time  value=$otime2  error=0.02
test tag=depth2 float=$depth value=$odepth2 error=0.15
test tag=rough2 float=$rough value=$orough2 error=0.30
