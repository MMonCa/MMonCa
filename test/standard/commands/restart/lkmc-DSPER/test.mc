
set otime1  11.1991
set odepth1 50.62     
set orough1 0.298317

set otime2  41.9355
set odepth2 45.983     
set orough2 1.90903

restart load=dump
report all

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

test tag=time1  float=$time  value=$otime1  error=0
test tag=depth1 float=$depth value=$odepth1 error=0
test tag=rough1 float=$rough value=$orough1 error=0

anneal temp=550 events=400000 time=5e3

set time  [extract time]
set depth [lindex [extract ac.mean] 0]
set rough [lindex [extract ac.stdev] 0]

test tag=time2  float=$time  value=$otime2  error=0.02
test tag=depth2 float=$depth value=$odepth2 error=0.02
test tag=rough2 float=$rough value=$orough2 error=0.02
