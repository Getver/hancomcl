Diversity_AB () {

    # METADATA=/work/meta_data.tsv

    # echo "[INFO] Running core diversity metrics"

    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny tree/rooted_tree.qza \
        --i-table filtered/denoise_noMT_table.qza \
        --p-sampling-depth 1109 \
        --m-metadata-file ${METADATA} \
        --output-dir diversity

    # echo "[INFO] Running alpha rarefaction"

    qiime diversity alpha-rarefaction \
        --i-table filtered/denoise_noMT_table.qza \
        --i-phylogeny tree/rooted_tree.qza \
        --p-max-depth 4000 \
        --m-metadata-file ${METADATA} \
        --o-visualization diversity/alpha_rarefaction.qzv

    # echo "[INFO] Exporting alpha rarefaction"

    qiime tools export \
        --input-path diversity/alpha_rarefaction.qzv \
        --output-path diversity/alpha

    # echo "[INFO] Exporting beta diversity PCoA results"

    mkdir -p diversity/txt

    for metric in bray_curtis jaccard unweighted_unifrac weighted_unifrac; do
        qiime tools export \
            --input-path diversity/${metric}_pcoa_results.qza \
            --output-path diversity/txt

        mv diversity/txt/ordination.txt \
           diversity/txt/${metric}_pcoa_results.txt
    done

    # echo "[INFO] Converting Emperor QZV to HTML"

    for qzv in diversity/*emperor*.qzv; do
        base=$(basename "$qzv" .qzv)
        mkdir -p "diversity/beta/${base}"

        tmpdir=$(mktemp -d)
        unzip -q -o "$qzv" -d "$tmpdir"

        uuid_dir=$(find "$tmpdir" -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -d "$uuid_dir/data" ]; then
            mv "$uuid_dir/data" "diversity/beta/${base}/data"
        fi

        rm -rf "$tmpdir"
    done

    # echo "[INFO] Diversity analysis completed"
}
