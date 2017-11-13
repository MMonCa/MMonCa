#Collect
echo "Collecting test information"
>test.tmp.results
for name in `find .`
do
	if [ -f $name/test.result ]; then
		cat $name/test.result >>test.tmp.results
	fi
done                                        
passed=`grep PASSED test.tmp.results | wc -l`
failed=`grep FAILED test.tmp.results | wc -l`
total=`expr $passed + $failed`
echo "$passed / $total tests PASSED" >>test.tmp.results
echo "$failed failure/s" >>test.tmp.results
cp test.tmp.results test.results
if [ $failed != 0 ]; then
	echo "--------- failed tests -------------" >>test.results
	grep FAILED test.tmp.results >>test.results
fi
rm test.tmp.results
cat test.results
echo "-- DONE --"
