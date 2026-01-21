Inspect_metadata () {

    qiime tools inspect-metadata "${MANIFEST}"

    if [[ "$SEQ_MODE" == "paired" ]]; then
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "${MANIFEST}" \
            --output-path "${OUTPUT_NAME}.qza" \
            --input-format PairedEndFastqManifestPhred33V2
    elif [[ "$SEQ_MODE" == "single" ]]; then
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "${MANIFEST}" \
            --output-path "${OUTPUT_NAME}.qza" \
            --input-format SingleEndFastqManifestPhred33V2
    fi

    qiime demux summarize \
        --i-data "${OUTPUT_NAME}.qza" \
        --o-visualization summary.qzv

    qiime tools export \
        --input-path summary.qzv \
        --output-path summary
}

    #############################################################################
    # 만약, fastq 파일에 barcode가 있다면 import 이렇게

    # mkdir -p demultiplex

    # ## demulti
    # qiime tools import \
    #   --type 'EMPPairedEndSequences' \
    #   --input-path ./pooled-fastq-directory/ \
    #   --output-path demultiplex/pooled_seqs.qza

    # qiime demux emp-paired \
    #   --i-seqs demultiplex/pooled_seqs.qza \
    #   --m-barcodes-file sample-metadata.tsv \
    #   --m-barcodes-column barcode-sequence \
    #   --o-per-sample-sequences demultiplex/${OUTPUT_NAME}_demulti.qza \
    #   --o-error-correction-details demultiplex/${OUTPUT_NAME}-details.qza

    # qiime demux summarize \
    #   --i-data demultiplex/${OUTPUT_NAME}_demulti.qza \
    #   --o-visualization demultiplex/${OUTPUT_NAME}_demulti_summary.qzv

    # qiime tools export \
    #   --input-path demultiplex/${OUTPUT_NAME}_demulti_summary.qzv \
    #   --output-path summary

#############################################################################



