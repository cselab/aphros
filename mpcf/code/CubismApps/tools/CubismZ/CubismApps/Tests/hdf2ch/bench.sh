TEST=$1
DATE=`date +%Y-%m-%d:%H:%M:%S`
LOG=/tmp/${DATE}.`echo ${TEST#./} | tr " " "-" `log.txt

for PARAM in 0.000001 0.000010 0.000100 0.001000 0.010000 0.100000 1.000000 
do
    echo "TEST with PARAM=$PARAM"
    echo "LOG AT $LOG"
    $TEST $PARAM
done | tee $LOG

cat $LOG | egrep "RES" | tee last-results.txt
#cat $LOG | egrep "RES" | awk '{ print $2, $4, $6, $8}' | tee last-results.txt

