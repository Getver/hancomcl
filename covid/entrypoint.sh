input_file=$1

day1_epoch=`date +%s`
day2_epoch=`stat -c %Y /files/corona_vari.txt`

T=$(($day1_epoch - $day2_epoch))

if [ ${T} -gt 2592000 ]; then
    python3 /opt/main/makeDB.py
fi

python3 /opt/main/admin.py
python3 /opt/main/run.py ${input_file}
