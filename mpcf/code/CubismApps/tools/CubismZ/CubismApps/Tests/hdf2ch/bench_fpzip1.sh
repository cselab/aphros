TEST=$1
DATE=`date +%Y-%m-%d:%H:%M:%S`
LOG=/tmp/${DATE}.`echo ${TEST#./} | tr " " "-" `log.txt
#RES=fpzip16_${DATE}.`echo ${TEST#./} | tr " " "-" `_res.txt
RES=`echo $2 | sed 's/\// /g' | awk '{print $NF}'`_fpzip1_res.txt

for PARAM in 31 30 28 24 22 20 18 16 14 12 10 #8
do
    echo "TEST with PARAM=$PARAM"
    echo "LOG AT $LOG"
    $TEST $2 $PARAM
done | tee $LOG


#$TEST $2 | tee $LOG

cat $LOG | egrep "RES" | tee $RES
#cat $LOG | egrep "RES" | tee last-results.txt
#cat $LOG | egrep "RES" | awk '{ print $2, $4, $6, $8}' | tee last-results.txt

