## QC : fastp
QC_fastp () {
    local SAMPLE="$1"
    local READ1="$2"
    local READ2="$3"

    mkdir -p fastp

    fastp \
        -i ${READ1} \
        -I ${READ2} \
        -o /home/data/fastp/${SAMPLE}_1.fastq.gz \
        -O /home/data/fastp/${SAMPLE}_2.fastq.gz \
        --n_base_limit 5 \
        --cut_front \
        --cut_tail \
        --average_qual 20 \
        --cut_mean_quality 20 \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --length_limit 400 \
        --detect_adapter_for_pe \
        --html fastp/${SAMPLE}_QC.html \
        --json fastp/${SAMPLE}_QC.json \
        --thread 10

        # --adapter_sequence <adapter_seq_R1> \
        # --adapter_sequence_r2 <adapter_seq_R2> \
}
# real    1m59.190s
# user    15m2.923s
# sys     0m8.570s