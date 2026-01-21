
qiime tools export \
  --input-path filtered/denoise_noMT_table.qza \
  --output-path filtered/tmp


/work/DB/picrust2-2.5.1/scripts/picrust2_pipeline.py \
  -s filtered/seq/denoise_noMT_seq.fasta \
  -i filtered/tmp/feature-table.biom \
  -o picrust/ \
  --processes 10 \
  --stratified \
  --per_sequence_contrib

