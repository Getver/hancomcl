Denoise_table () {

    mkdir -p filtered

    ## Run
    # noise 제거(전처리), Run, TEST3, TEST4
    # qiime dada2 denoise-paired \
    #     --i-demultiplexed-seqs ${OUTPUT_NAME}.qza \
    #     --p-trim-left-f 0 \
    #     --p-trim-left-r 0 \
    #     --p-trunc-len-f 280 \
    #     --p-trunc-len-r 240 \
    #     --p-max-ee-f 5 \
    #     --p-max-ee-r 5 \
    #     --p-trunc-q 3 \
    #     --o-table filtered/denoise_table.qza \
    #     --o-representative-sequences filtered/denoise_seq.qza \
    #     --o-denoising-stats filtered/denoise_stat.qza \
    #     --p-chimera-method consensus \
    #     --p-n-threads 50 \
    #     --verbose
    
    ## TEST1, TEST2
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${OUTPUT_NAME}.qza \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-trunc-len-f 145 \
        --p-trunc-len-r 140 \
        --p-max-ee-f 5 \
        --p-max-ee-r 5 \
        --p-trunc-q 2 \
        --o-table filtered/denoise_table.qza \
        --o-representative-sequences filtered/denoise_seq.qza \
        --o-denoising-stats filtered/denoise_stat.qza \
        --p-chimera-method consensus \
        --p-n-threads 20 \
        --verbose

    qiime metadata tabulate \
    --m-input-file filtered/denoise_stat.qza \
    --o-visualization filtered/denoise_stat.qzv

    ## stat 파일로 변환
    qiime tools export \
    --input-path filtered/denoise_stat.qza \
    --output-path filtered/txt/

    mv filtered/txt/stats.tsv filtered/txt/denoise_stat.tsv

    # table 파일로 변환
    qiime tools export \
    --input-path filtered/denoise_table.qza \
    --output-path filtered/

    qiime feature-table summarize \
    --i-table filtered/denoise_table.qza \
    --o-visualization filtered/denoise_table.qzv \
    --m-sample-metadata-file ${MANIFEST}

    biom convert \
        -i filtered/feature-table.biom \
        -o filtered/txt/denoise_table.tsv \
        --to-tsv
    
    rm filtered/feature-table.biom

    qiime tools export \
    --input-path filtered/denoise_table.qzv \
    --output-path filtered/denoiseTable

    qiime tools export \
    --input-path filtered/denoise_seq.qza \
    --output-path filtered/seq

    mv filtered/seq/dna-sequences.fasta filtered/seq/denoise_seq.fasta
}
# real    68m50.637s
# user    550m50.924s
# sys     3m25.631s