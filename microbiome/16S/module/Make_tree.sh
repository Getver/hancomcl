Make_tree () {

    mkdir -p tree


    ## 계통수 구축하기
    qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences filtered/denoise_noMT_seq.qza \
    --o-alignment tree/aligned.qza \
    --o-masked-alignment tree/masked.qza \
    --o-tree tree/unrooted.qza \
    --o-rooted-tree tree/rooted_tree.qza

    ## tree Newick export
    qiime tools export \
    --input-path tree/rooted_tree.qza \
    --output-path tree

    # ----------------------------------------------
    # 2. Empress 시각화 (Visualization)
    # ----------------------------------------------

  ## 속 수준
  # qiime taxa collapse \
  #   --i-table filtered/denoise_noMT_table.qza \
  #   --i-taxonomy taxonomy/taxonomy.qza \
  #   --p-level 6 \
  #   --o-collapsed-table filtered/table_L6.qza



    # qiime empress community-plot \
    # --i-tree tree/rooted_tree.qza \
    # --i-feature-table filtered/table_L6.qza \
    # --m-sample-metadata-file ${MANIFEST} \
    # --m-feature-metadata-file taxonomy/taxonomy.qza \
    # --p-filter-missing-features \
    # --o-visualization tree/empress_tree.qzv

    # ----------------------------------------------
    # 3. HTML 변환 (웹 리포트용)
    # ----------------------------------------------
    
    # qiime tools export \
    # --input-path tree/empress_tree.qzv \
    # --output-path tree/empress_html


    # qiime phylogeny view-tree \
    # --i-tree tree/tree.qza \
    # --o-visualization tree/tree.qzv
    # qiime empress tree-plot \
    # --i-tree tree/tree.qza \
    # --i-feature-table filtered/denoise_noMT.qza \
    # --i-taxonomy taxonomy/taxonomy.qza \
    # --o-visualization tree/tree-empress.qzv
    }



