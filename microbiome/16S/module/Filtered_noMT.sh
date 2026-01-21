Filtered_noMT () {
    # local INPUT_FILE="$1"
    
    ## 미토콘드리아 없애기(table)
    qiime taxa filter-table \
    --i-table filtered/denoise_table.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --p-mode contains \
    --p-include p__ \
    --p-exclude 'p__;,Chloroplast,Mitochondria' \
    --o-filtered-table filtered/denoise_noMT_table.qza

    ## 미토콘드리아 없애기(seq)
    qiime taxa filter-seqs  \
    --i-sequences filtered/denoise_seq.qza\
    --i-taxonomy taxonomy/taxonomy.qza \
    --p-mode contains \
    --p-include p__ \
    --p-exclude 'p__;,Chloroplast,Mitochondria' \
    --o-filtered-sequences  filtered/denoise_noMT_seq.qza

    ## 이게 샘플 QC인데
    qiime feature-table filter-samples \
    --i-table filtered/denoise_noMT_table.qza \
    --p-min-frequency 1000 \
    --o-filtered-table filtered/denoise_noMT_table.qza

    qiime feature-table summarize \
        --i-table filtered/denoise_noMT_table.qza \
        --o-visualization filtered/denoise_noMT_table.qzv \
        --m-sample-metadata-file ${MANIFEST}

    # table 파일 생성(biom)
    qiime tools export \
        --input-path filtered/denoise_noMT_table.qza \
        --output-path filtered

    biom convert \
        -i filtered/feature-table.biom \
        -o filtered/txt/denoise_noMT_table.tsv \
        --to-tsv

    rm filtered/feature-table.biom
    
    qiime tools export \
    --input-path filtered/denoise_noMT_table.qzv \
    --output-path filtered/noMTable

    qiime tools export \
    --input-path filtered/denoise_noMT_seq.qza \
    --output-path filtered/seq

    mv filtered/seq/dna-sequences.fasta filtered/seq/denoise_noMT_seq.fasta
}
# real    0m40.795s
# user    0m51.972s
# sys     0m1.791s