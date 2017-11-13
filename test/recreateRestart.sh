#!/bin/sh
# 16 means a maximum of 8 CPUs busy (each one is counted as two because of the time command)
if [ "A$1" = "A" ]; then
   version=Obj_g++
else
   version=$1
fi
echo "Running for configuration $version"
current=`pwd`
export MCPATH=$current/../config
echo "using $MCPATH"
for name in `find .`
do
	if [ -d $name -a -f $name/generate.mc ]; then
		cd $name
		$current/../$version/mmonca generate.mc &
		cd $current
	fi
done
