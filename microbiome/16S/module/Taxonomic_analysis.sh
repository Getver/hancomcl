Taxonomic_analysis () {

    mkdir -p taxonomy

    ## Taxonomic analysis
    # Debug info has been saved to /tmp/qiime2-q2cli-err-u41dmn7z.log
    # --i-classifier ~/references/microbiome/silva-138-99-nb-weighted-classifier.qza \
    qiime feature-classifier classify-sklearn \
    --i-classifier /work/DB/silva-138.2-ssu-nr99-classifier.qza \
    --i-reads filtered/denoise_seq.qza \
    --o-classification taxonomy/taxonomy.qza

    ## tsv 생성
    qiime tools export \
    --input-path taxonomy/taxonomy.qza \
    --output-path taxonomy

    mv taxonomy/taxonomy.tsv taxonomy/taxonomy_silva.tsv
    
}

# real    1m2.338s
# user    1m16.408s
# sys     0m11.594s