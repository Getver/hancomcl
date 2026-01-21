Taxonomy_kraken () {
    local SAMPLE="$1"

    mkdir -p kraken

    kraken2 \
        --db /work/DB/kraken2_db \
        --threads 10 \
        --report kraken/${SAMPLE}.report \
        --output kraken/${SAMPLE}.kraken \
        --paired host_removed/${SAMPLE}_host_removed_R1.fastq.gz host_removed/${SAMPLE}_host_removed_R2.fastq.gz

    # kraken2 \
    #     --db /work/DB/kraken2_db \
    #     --threads 10 \
    #     --report kraken/${SAMPLE}_mpa.report \
    #     --use-mpa-style \
    #     --output kraken/${SAMPLE}_2.kraken \
    #     --paired host_removed/${SAMPLE}_host_removed_R1.fastq.gz host_removed/${SAMPLE}_host_removed_R2.fastq.gz
}
# real    6m10.156s
# user    3m59.606s
# sys     1m10.041s