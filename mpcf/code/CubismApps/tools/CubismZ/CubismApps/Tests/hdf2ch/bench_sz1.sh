TEST=$1
DATE=`date +%Y-%m-%d:%H:%M:%S`
LOG=/tmp/${DATE}.`echo ${TEST#./} | tr " " "-" `log.txt
#RES=sz_${DATE}.`echo ${TEST#./} | tr " " "-" `_res.txt
RES=`echo $2 | sed 's/\// /g' | awk '{print $NF}'`_sz1_res.txt

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-1.8.0/jre/lib/amd64/server

#for PARAM in 0.000000 0.000001 0.000010 0.000100 0.001000 0.010000 0.100000 1 10 50 100 200 500 1000
for PARAM in 0 0.00000001 0.0000001 0.000001 0.000010 0.000100 0.001000 0.010000 0.100000 1 10 50 100 200 500 1000
do
    echo "TEST with PARAM=$PARAM"
    echo "LOG AT $LOG"
    $TEST $2 $PARAM
done | tee $LOG

cat $LOG | egrep "RES" | tee $RES
#cat $LOG | egrep "RES" | tee last-results.txt
#cat $LOG | egrep "RES" | awk '{ print $2, $4, $6, $8}' | tee last-results.txt

