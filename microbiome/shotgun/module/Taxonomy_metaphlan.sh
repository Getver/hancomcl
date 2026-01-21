## Taxonomy profiling : metaphlan
Taxonomy_metaphlan () {
    local SAMPLE="$1"

    mkdir -p metaphlan

    metaphlan \
        host_removed/${SAMPLE}_merged.fastq.gz \
        --input_type fastq \
        --nproc 10 \
        -o metaphlan/${SAMPLE}_profile.txt
        # --bowtie2out metaphlan/${SAMPLE}.bowtie2.bz2

    # metaphlan \
    #     host_removed/${SAMPLE}_merged.fastq.gz.bowtie2out.txt \
    #     --input_type bowtie2out \
    #     --nproc 10 \
    #     -o metaphlan/${SAMPLE}_profile.txt
}

# real    2m34.900s
# user    21m8.361s
# sys     0m25.825s