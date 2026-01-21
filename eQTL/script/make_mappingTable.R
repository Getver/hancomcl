# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=38)

genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)

mapping_df <- genes %>% select(ensembl_gene_id, hgnc_symbol, chromosome_name, start_position, end_position) %>% unique()

write.table(mapping_df, '../text/mapping_table.txt', sep='\t', quote=FALSE, row.names=FALSE)

