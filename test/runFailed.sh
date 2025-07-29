#!/bin/sh
# 16 means a maximum of 8 CPUs busy (each one is counted as two because of the time command)
if [ "A$1" = "A" ]; then
   version=Obj_Test
else
   version=$1
fi
echo "Running only FAILED ones for configuration $version"
current=`pwd`
export MCPATH=$current/../config
echo "using $MCPATH"
rm -f tests.started
for name in `find .`
do
	if [ -d $name -a -f $name/test.mc ]; then
		failed=`grep FAILED $name/test.result | wc -l`
		if [ $failed != 0 ]; then
			./runSingle.sh $current $name $version &
		else
			cat $name/test.result
		fi
	fi
	while [ `ps aux | grep mmonca | grep test | wc -l` -ge 2 ]; do
		sleep 1
	done
done
#see if everyone has finished
while [ `ps aux | grep mmonca | grep test | wc -l` -ne 0 ]; do
	sleep 1
done
./collect.sh
