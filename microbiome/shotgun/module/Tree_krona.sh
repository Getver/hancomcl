Tree_krona () {
    local SAMPLE="$1"

    ktImportTaxonomy \
    -t 5 \
    -m 3 \
    -o kraken/${SAMPLE}.krona.html \
    kraken/${SAMPLE}.report


    # kraken2 \
    #     --db /work/DB/kraken2_db \
    #     --threads 10 \
    #     --report kraken/${SAMPLE}_mpa.report \
    #     --use-mpa-style \
    #     --output kraken/${SAMPLE}_2.kraken \
    #     --paired host_removed/${SAMPLE}_host_removed_R1.fastq.gz host_removed/${SAMPLE}_host_removed_R2.fastq.gz
}