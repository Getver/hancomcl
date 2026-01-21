Denoise_table () {

    mkdir -p filtered



    if [[ "$SEQ_MODE" == "paired" ]]; then
        qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${OUTPUT_NAME}.qza \
        --p-front-f CCTACGGGNGGCWGCAG \
        --p-front-r GACTACHVGGGTATCTAATCC \
        --p-discard-untrimmed \
        --o-trimmed-sequences demux_noprimer.qza


        # qiime dada2 denoise-paired \
        #     --i-demultiplexed-seqs ${OUTPUT_NAME}.qza \
        #     --p-trim-left-f 0 \
        #     --p-trim-left-r 0 \
        #     --p-trunc-len-f 145 \
        #     --p-trunc-len-r 140 \
        #     --p-max-ee-f 5 \
        #     --p-max-ee-r 5 \
        #     --p-trunc-q 2 \
        #     --o-table filtered/denoise_table.qza \
        #     --o-representative-sequences filtered/denoise_seq.qza \
        #     --o-denoising-stats filtered/denoise_stat.qza \
        #     --p-chimera-method consensus \
        #     --p-n-threads 20 \
        #     --verbose

        ## Run
        # noise 제거(전처리), Run, TEST3, TEST4
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs demux_noprimer.qza \
            --p-trim-left-f 0 \
            --p-trim-left-r 0 \
            --p-trunc-len-f 280 \
            --p-trunc-len-r 240 \
            --p-max-ee-f 5 \
            --p-max-ee-r 5 \
            --p-trunc-q 3 \
            --o-table filtered/denoise_table.qza \
            --o-representative-sequences filtered/denoise_seq.qza \
            --o-denoising-stats filtered/denoise_stat.qza \
            --p-chimera-method consensus \
            --p-n-threads 12 \
            --verbose
        
        qiime tools export \
        --input-path filtered/denoise_stat.qza \
        --output-path filtered/txt/
        
        qiime metadata tabulate \
        --m-input-file filtered/denoise_stat.qza \
        --o-visualization filtered/denoise_stat.qzv


        # qiime demux summarize \
        # --i-data filtered/denoise_stat.qza \
        # --o-visualization filtered/denoise_stat.qzv

    elif [[ "$SEQ_MODE" == "single" ]]; then
        # =========================
        # Nanopore (VSEARCH OTU)
        # =========================
        qiime cutadapt trim-single \
        --i-demultiplexed-sequences ${OUTPUT_NAME}.qza \
        --p-minimum-length 1300 \
        --o-trimmed-sequences filtered/denoise_seq_tmp.qza

        qiime vsearch dereplicate-sequences \
            --i-sequences filtered/denoise_seq_tmp.qza \
            --o-dereplicated-table filtered/denoise_table_tmp.qza \
            --o-dereplicated-sequences filtered/rep_seqs.qza

        qiime vsearch cluster-features-de-novo \
            --i-table filtered/denoise_table_tmp.qza \
            --i-sequences filtered/rep_seqs.qza \
            --p-perc-identity 0.95 \
            --p-threads 50 \
            --o-clustered-table filtered/denoise_table.qza \
            --o-clustered-sequences filtered/denoise_seq.qza
    fi

    # ---------------------------
    # stats 요약
    # ---------------------------

    # qiime demux summarize \
    #     --i-data filtered/denoise_seq_tmp.qza \
    #     --o-visualization filtered/denoise_seq_tmp.qzv

    # qiime metadata tabulate \
    #     --m-input-file filtered/denoise_seq.qza \
    #     --o-visualization filtered/denoise_seq.qzv


    # qiime tools export \
    #     --input-path filtered/denoise_stat.qza \
    #     --output-path filtered/txt/

    # mv filtered/txt/stats.csv filtered/txt/denoise_seq.tsv

    # ---------------------------
    # table export
    # ---------------------------
    # qiime feature-table filter-features \
    # --i-table filtered/denoise_table.qza \
    # --p-min-frequency 10 \
    # --o-filtered-table filtered/denoise_table.qza

    qiime tools export \
        --input-path filtered/denoise_table.qza \
        --output-path filtered/

    qiime feature-table summarize \
        --i-table filtered/denoise_table.qza \
        --o-visualization filtered/denoise_table.qzv \
        --m-sample-metadata-file "${MANIFEST}"
    
    # mkdir -p filtered/txt

    biom convert \
        -i filtered/feature-table.biom \
        -o filtered/txt/denoise_table.tsv \
        --to-tsv

    # rm filtered/feature-table.biom

    qiime tools export \
        --input-path filtered/denoise_table.qzv \
        --output-path filtered/denoiseTable

    # ---------------------------
    # representative sequences
    # ---------------------------
    qiime tools export \
        --input-path filtered/denoise_seq.qza \
        --output-path filtered/seq

    mv filtered/seq/dna-sequences.fasta filtered/seq/denoise_seq.fasta

}
