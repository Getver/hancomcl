## FASTQ 이름 변경 및 병합
Preprocessing_FASTQ () {
    local SAMPLE="$1"

    ## 나중을 생각하여 넣는 것이 좋을 것 같은 코드임
    mv host_removed/${SAMPLE}_host_removed.1 host_removed/${SAMPLE}_host_removed_R1.fastq.gz
    mv host_removed/${SAMPLE}_host_removed.2 host_removed/${SAMPLE}_host_removed_R2.fastq.gz

    ## 얘는 꼭 필요
    zcat host_removed/${SAMPLE}_host_removed_R1.fastq.gz host_removed/${SAMPLE}_host_removed_R2.fastq.gz | gzip > host_removed/${SAMPLE}_merged.fastq.gz
}
# real    40m29.538s
# user    44m52.011s
# sys     0m19.707s