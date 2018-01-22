TEST=$1
DATE=`date +%Y-%m-%d:%H:%M:%S`
LOG=/tmp/${DATE}.`echo ${TEST#./} | tr " " "-" `log.txt
#RES=fpzip24_${DATE}.`echo ${TEST#./} | tr " " "-" `_res.txt
RES=`echo $2 | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt

$TEST $2 | tee $LOG

cat $LOG | egrep "RES" | tee $RES
#cat $LOG | egrep "RES" | tee last-results.txt
#cat $LOG | egrep "RES" | awk '{ print $2, $4, $6, $8}' | tee last-results.txt

