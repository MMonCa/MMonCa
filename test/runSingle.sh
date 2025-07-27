#!/bin/sh
# args: current directory, directory of test, binary.
# example
#  ./runSingle.sh `pwd` okmc/defects/IVCluster/Rec Release
if [ "A$1" = "A" ]; then
	echo "I need the current directory "
fi
if [ "A$2" = "A" ]; then
	echo "I need the directory where the test is"
fi
cd $2
if [ -f errors ]; then
	rm errors
fi
echo "$2/test.mc" >>$1/tests.started
if time -o time -f "%U" $1/../$3/mmonca test.mc >/dev/null 2> test.errors; then
	echo "$2/test.mc `cat time` PASSED" >test.result
	echo "$2/test.mc `cat time` PASSED"
	rm test.errors
else
	echo "$2/test.mc FAILED" >test.result
	echo "$2/test.mc FAILED"
	tail test.errors
	mv test.errors errors
fi
cd $1
