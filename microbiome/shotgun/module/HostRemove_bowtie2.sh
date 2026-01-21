## Host remove : bowtie2
HostRemove_bowtie2 () {
    local SAMPLE="$1"

    mkdir -p host_removed

    bowtie2 \
        --threads 10 \
        -x /work/DB/host_genome/hg38 \
        -1 /home/data/fastp/${SAMPLE}_1.fastq.gz \
        -2 /home/data/fastp/${SAMPLE}_2.fastq.gz \
        --very-sensitive-local \
        --un-conc-gz host_removed/${SAMPLE}_host_removed > host_removed/${SAMPLE}.sam \
        2> host_removed/${SAMPLE}_bowtie2.log
}
# real    23m43.007s
# user    143m40.971s
# sys     65m25.628s