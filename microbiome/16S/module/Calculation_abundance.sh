Calculation_abundance () {

  mkdir -p abundance

  # # 1. Feature table를 genus 수준으로 합치기
  # qiime taxa collapse \
  #   --i-table filtered/${OUTPUT_NAME}_denoise_noMT_table.qza \
  #   --i-taxonomy taxonomy/${OUTPUT_NAME}_taxonomy.qza \
  #   --p-level 6 \
  #   --o-collapsed-table taxonomy/table-genus.qza

  # # 2. Genus 테이블에서 상대 풍부도 계산
  # qiime feature-table relative-frequency \
  #   --i-table taxonomy/table-genus.qza \
  #   --o-relative-frequency-table abundance/${OUTPUT_NAME}_genus_abundance.qza

  # # 3. TSV로 export
  # qiime tools export \
  #   --input-path abundance/${OUTPUT_NAME}_genus_abundance.qza \
  #   --output-path abundance

  # # 4. biom → tsv 변환
  # biom convert -i abundance/feature-table.biom \
  #             -o abundance/${OUTPUT_NAME}_genus_abundance.tsv \
  #             --to-tsv


  qiime feature-table relative-frequency \
    --i-table filtered/denoise_noMT_table.qza \
    --o-relative-frequency-table abundance/abundance.qza

  qiime tools export \
    --input-path abundance/abundance.qza \
    --output-path abundance

  biom convert -i abundance/feature-table.biom \
              -o abundance/abundance.tsv \
              --to-tsv
}