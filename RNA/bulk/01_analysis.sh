## STAR v2 사용하여 STAR Alignment
# https://academic.oup.com/bioinformatics/article/29/1/15/272537?login=false

# docker run -it -v /disk0/sm/bulk/kobic:/root -v /disk0/references/hg38/UCSC/RNA_for_kobic/gencode:/data test:0.1



DATA=`cat /disk0/sm/bulk/test/fastq.txt`

echo ${DATA}

for DA in ${DATA}
do

    FILE=`echo ${DA} | sed 's/,/ /g'`

    D=`echo ${FILE} | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}' | awk -F '_' '{FS="_"; OFS="_"} {$NF=""; print $0}' | rev | cut -c2- | rev`
    echo ${FILE} ${D}

    # fastqc -o /disk0/sm/bulk/test/fastqc ${FILE}
    
    # mkdir /disk0/sm/bulk/test/output/${D}

    STAR \
    --readFilesIn ${FILE} \
    --outSAMattrRGline "ID:${D} SM:${D}" \
    --genomeDir /disk0/references/hg38/UCSC/RNA/all_chromosomes/gencode \
    --readFilesCommand zcat \
    --runThreadN 99 \
    --twopassMode Basic \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterType BySJout \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 \
    --limitSjdbInsertNsj 1200000 \
    --outFileNamePrefix /disk0/sm/bulk/test/output/${D}/${D}_ \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs None \
    --alignSoftClipAtReferenceEnds Yes \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM SortedByCoordinate Unsorted \
    --outSAMunmapped Within \
    --genomeLoad NoSharedMemory \
    --chimSegmentMin 15 \
    --chimJunctionOverhangMin 15 \
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
    --chimOutJunctionFormat 1 \
    --chimMainSegmentMultNmax 1 \
    --outSAMattributes NH HI AS nM NM ch

    rsem-calculate-expression \
    --no-bam-output \
    -p 18 \
    --forward-prob 0 \
    --estimate-rspd \
    --alignments \
    --paired-end \
    /disk0/sm/bulk/test/output/${D}/${D}_Aligned.toTranscriptome.out.bam \
    /disk0/references/hg38/UCSC/RNA/all_chromosomes/gencode/rsem/hg38_rsem \
    /disk0/sm/bulk/test/output/${D}/${D}_
    
done



