# ONESAMPLE_PIPELINE () {
#     local SAMPLE="$1"
#     local SAMPLE1="$2"
#     local SAMPLE2="$3"

#     ## QC : fastp
#     QC_fastp ${SAMPLE} ${SAMPLE1} ${SAMPLE2}

#     ## Host remove : bowtie2
#     HostRemove_bowtie2 ${SAMPLE}
#     Preprocessing_FASTQ ${SAMPLE}

#     ## Taxonomy profiling : metaphlan
#     Taxonomy_metaphlan ${SAMPLE}

#     ## Functional profiling : HUMAnN3
#     Functional_humann ${SAMPLE}

#     ## Taxonomy profiling : kraken
#     Taxonomy_kraken ${SAMPLE}

#     ## Abundance profiling : bracken
#     Abundance_bracken ${SAMPLE}
# }

ONESAMPLE_PIPELINE () {
    local SAMPLE="$1"
    local SAMPLE1="$2"
    local SAMPLE2="$3"
    local LOGFILE="/home/data/logs/${SAMPLE}_runtime.log"

    mkdir -p logs
    echo "=== Runtime log for ${SAMPLE} ===" > "$LOGFILE"

    log_runtime() {
        local step_name="$1"
        shift
        local start_time=$(date +%s)
        "$@"   # 실행
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        printf "[%s] %02d:%02d:%02d\n" \
            "$step_name" $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60)) \
            >> "$LOGFILE"
    }

    log_runtime "QC_fastp" QC_fastp "$SAMPLE" "$SAMPLE1" "$SAMPLE2"
    log_runtime "HostRemove_bowtie2" HostRemove_bowtie2 "$SAMPLE"
    log_runtime "Preprocessing_FASTQ" Preprocessing_FASTQ "$SAMPLE"
    log_runtime "Taxonomy_metaphlan" Taxonomy_metaphlan "$SAMPLE"
    log_runtime "Functional_humann" Functional_humann "$SAMPLE"
    log_runtime "Taxonomy_kraken" Taxonomy_kraken "$SAMPLE"
    log_runtime "Abundance_bracken" Abundance_bracken "$SAMPLE"
    log_runtime "Tree_krona" Tree_krona "$SAMPLE"


    echo "=== Pipeline finished for ${SAMPLE} ===" >> "$LOGFILE"
}

