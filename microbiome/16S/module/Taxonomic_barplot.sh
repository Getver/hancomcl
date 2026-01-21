Taxonomic_barplot () {
    # local INPUT_FILE="$1"

    qiime tools import \
    --input-path taxonomy/ncbi_qiime2.tsv \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --output-path taxonomy/taxonomy_blast.qza
    
    ## tsv 생성
    qiime tools export \
    --input-path taxonomy/taxonomy_blast.qza \
    --output-path taxonomy

    mv taxonomy/taxonomy.tsv taxonomy/taxonomy_blast.tsv

    # ----- silva랑 비교 후 

    ## 메타데이터(meta.txt) 없이 하려면 그냥 file_list.txt 써도 됨
    qiime taxa barplot \
    --i-table filtered/denoise_noMT_table.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --m-metadata-file ${MANIFEST} \
    --o-visualization taxonomy/taxonomyBar.qzv
                # --m-metadata-file ${MANIFEST} \

    qiime tools export \
        --input-path taxonomy/taxonomyBar.qzv \
        --output-path taxonomy/taxonomyBar
}


    # ## 메타데이터(meta.txt) 없이 하려면 그냥 file_list.txt 써도 됨
    # qiime taxa barplot \
    #     --i-table filtered/denoise_noMT_table.qza \
    #     --i-taxonomy taxonomy/taxonomy.qza \
    #     --m-metadata-file ${MANIFEST} \
    #     --o-visualization taxonomy/taxonomyBar.qzv

    # qiime tools export \
    #     --input-path taxonomy/taxonomyBar.qzv \
    #     --output-path taxonomy/taxonomyBar