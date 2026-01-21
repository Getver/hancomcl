Abundance_bracken () {
    local SAMPLE="$1"

    bracken \
        -d /work/DB/kraken2_db \
        -i kraken/${SAMPLE}.report \
        -o kraken/${SAMPLE}_bracken.txt \
        -r 150 \
        -l S
}