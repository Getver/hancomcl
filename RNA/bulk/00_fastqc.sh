Folder='/disk0/sm/bulk/2024/1114'
Input='/disk2/hancomcl/mhchoi/241113_DNAlink_FASTQ'
#  multiqc fastqc/*_fastqc.zip

# 가천대
/disk0/sm/WES/tool/FastQC/fastqc -o ${Folder}/fastqc -t 10 ${Input}/*.fastq.gz


# https://www.incodom.kr/QC

html=`ls ${Folder}/fastqc/*.html`

for i in $html
do
    # echo ${i}
    name=`ls $i | awk -F '/' '{print $NF}' | awk -F '_' '{print $1"_"$2}'`
    # echo ${name}
    read=`cat $i | sed -e 's/<tr>/\n/g' | grep 'Total Sequences' | awk -F '[><]' '{print $7}'`
    length=`cat $i | sed -e 's/<tr>/\n/g' | grep 'Sequence length</td>' | awk -F '[><]' '{print $7}'`
    # echo ${read} ${length}
    read_length=`expr ${read} \* ${length} \* 2`
    list=`echo $name $read_length | awk -F ' ' '{print $1"\t"$2}' | sort | uniq`
    GC=`cat $i | sed -e 's/<tr>/\n/g' | grep '%GC</td>' | awk -F '[><]' '{print $7}'`
    echo $list  ${length} ${read} $GC
done > ${Folder}/total_sequence.txt

# (fastqc file unzip)
# for file in `ls ${Folder}/fastqc/*.zip`; do unzip "${file}"; done

DATA=`ls ${Folder}/fastqc/*.zip | awk -F '.' '{print $1}' | awk -F '/' '{print $NF}' | awk -F '_' '{print $1"_"$2}' | sort | uniq`

for D in ${DATA}
do
    file1=`echo ${Folder}/fastqc/${D}_R1_001_fastqc/fastqc_data.txt`
    file2=`echo ${Folder}/fastqc/${D}_R2_001_fastqc/fastqc_data.txt`
    total1=`grep -A36 'Per sequence quality scores' ${file1} | sed 1,2d | awk '{ sum += $2 } END { print sum }'`
    Q201=`grep -A36 'Per sequence quality scores' ${file1} | sed 1,2d | awk '{ if ($1 > 19) {sum += $2} } END { print sum }'` 
    Q301=`grep -A36 'Per sequence quality scores' ${file1} | sed 1,2d | awk '{ if ($1 > 29) {sum += $2} } END { print sum }'`
    total2=`grep -A36 'Per sequence quality scores' ${file2} | sed 1,2d | awk '{ sum += $2 } END { print sum }'`
    Q202=`grep -A36 'Per sequence quality scores' ${file2} | sed 1,2d | awk '{ if ($1 > 19) {sum += $2} } END { print sum }'` 
    Q302=`grep -A36 'Per sequence quality scores' ${file2} | sed 1,2d | awk '{ if ($1 > 29) {sum += $2} } END { print sum }'`
    dream=`echo ${D} ${Q201} ${Q301} ${total1} ${Q202} ${Q302} ${total2} | awk -F ' ' '{print $1,$2/$4*100,$3/$4*100,$5/$7*100,$6/$7*100}'`
    echo ${dream}

done > ${Folder}/total_qulity.txt

python3 fastqc.py ${Folder}