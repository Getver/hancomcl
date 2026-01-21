# https://pmc.ncbi.nlm.nih.gov/articles/PMC7365738/

# mkdir -p picrust



picrust2_pipeline.py \
  -s filtered/seq/denoise_noMT_seq.fasta \
  -i filtered/feature-table.biom \
  -o picrust/ \
  -p 16 \
  -t epa-ng \
  --no_regroup \
  --stratified \
  --per_sequence_contrib \
  --verbose





# # === setup.R ===
# required <- c("edgeR","limma","pheatmap","RColorBrewer","ggplot2","data.table","tidyr","dplyr", "R.utils", "pheatmap")
# install_if_missing <- function(pkgs){
#   to_inst <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
#   if(length(to_inst)) install.packages(to_inst, repos="https://cloud.r-project.org")
# }
# install_if_missing(required)

# library(edgeR); library(limma); library(pheatmap)
# library(RColorBrewer); library(ggplot2)
# library(data.table); library(tidyr); library(dplyr)

# # --- 파일 경로 (환경에 맞게 수정) ---
# ko_table_fp      <- "picrust/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"        # 또는 gz 제거/경로 확인
# path_table_fp    <- "picrust/pathways_out/path_abun_unstrat.tsv.gz"
# metadata_fp      <- "/work/meta_data.tsv"   # 너의 메타데이터 파일 경로

# # 읽기 helper (gz or plain)
# fread_maybe_gz <- function(fp){
#   if(file.exists(fp)) return(data.table::fread(fp))
#   if(file.exists(paste0(fp,".gz"))) return(data.table::fread(paste0(fp,".gz")))
#   stop("파일을 찾을 수 없습니다: ", fp)
# }


# ko <- fread_maybe_gz(ko_table_fp)   # 첫 column KO id, 다음이 샘플들
# # PICRUSt2 metagenome_unstrat 형식: 첫 col = feature (KO), 다음이 sample columns
# # 행이 KO, 열이 samples
# # 변환: data.table -> matrix (rownames=KO)
# ko_mat <- as.data.frame(ko)
# rownames(ko_mat) <- ko_mat[[1]]; ko_mat[[1]] <- NULL
# ko_mat <- as.matrix(ko_mat)

# meta <- fread_maybe_gz(metadata_fp)
# meta <- as.data.frame(meta)
# rownames(meta) <- meta[,1]   # assume first col is SampleID

# # Match sample order
# samps <- intersect(colnames(ko_mat), rownames(meta))
# ko_mat <- ko_mat[, samps, drop=FALSE]
# meta <- meta[samps, , drop=FALSE]

# # filter low-abundance KOs and choose top variable
# keep <- rowSums(ko_mat) > 1e-6    # remove zero rows
# ko_mat_f <- ko_mat[keep, ,drop=FALSE]
# varKO <- apply(ko_mat_f, 1, var)
# topN <- 100
# topKOs <- names(sort(varKO, decreasing=TRUE))[1:topN]
# mat_top <- ko_mat_f[topKOs, ]

# # scale rows for heatmap
# mat_z <- t(scale(t(log10(mat_top + 1))))

# # annotation for columns (example: Group column in metadata)
# anno_col <- data.frame(Group = meta$Group)
# rownames(anno_col) <- rownames(meta)

# pheatmap(mat_z,
#          cluster_rows=TRUE, cluster_cols=TRUE,
#          annotation_col = anno_col,
#          show_rownames = TRUE, show_colnames = TRUE,
#          fontsize_row = 6,
#          main = paste0("Top ", topN, " variable KOs (log10+1, row-z)"))


# # === pca_KO.R ===
# source("setup.R")

# ko <- fread_maybe_gz(ko_table_fp)
# ko_df <- as.data.frame(ko); rownames(ko_df) <- ko_df[[1]]; ko_df[[1]] <- NULL
# ko_mat <- as.matrix(ko_df)
# meta <- fread_maybe_gz(metadata_fp); meta <- as.data.frame(meta); rownames(meta) <- meta[,1]

# samps <- intersect(colnames(ko_mat), rownames(meta))
# ko_mat <- ko_mat[, samps]
# meta <- meta[samps, , drop=FALSE]

# # Normalize: TMM (edgeR) -> logCPM
# dge <- DGEList(counts = round(ko_mat))   # PICRUSt predicted abundances: approximate counts -> round
# dge <- calcNormFactors(dge, method="TMM")
# lcpm <- cpm(dge, log=TRUE, prior.count=1)

# # PCA
# pca <- prcomp(t(lcpm), center=TRUE, scale.=FALSE)
# pcvar <- round(100 * summary(pca)$importance[2, 1:2], 1)

# pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], meta)
# library(ggplot2)
# ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) +
#   geom_point(size=3) +
#   labs(x=paste0("PC1 (",pcvar[1], "%)"), y=paste0("PC2 (",pcvar[2], "%)"),
#        title="PCA on KO logCPM") +
#   theme_minimal()



# # === pathway_enrichment.R ===
# source("setup.R")

# path <- fread_maybe_gz(path_table_fp)
# # path table: first col = pathway id/name, next cols = samples
# path_df <- as.data.frame(path); rownames(path_df) <- path_df[[1]]; path_df[[1]] <- NULL
# path_mat <- as.matrix(path_df)

# meta <- fread_maybe_gz(metadata_fp); meta <- as.data.frame(meta); rownames(meta) <- meta[,1]
# samps <- intersect(colnames(path_mat), rownames(meta))
# path_mat <- path_mat[, samps]
# meta <- meta[samps, , drop=FALSE]

# # Design: 예, Group (two-group 비교 가정: Group has two levels A/B)
# meta$Group <- factor(meta$Group)
# design <- model.matrix(~ 0 + meta$Group)
# colnames(design) <- levels(meta$Group)

# # limma-voom
# dge <- DGEList(counts = round(path_mat))
# dge <- calcNormFactors(dge)
# v <- voom(dge, design=design, plot=FALSE)
# fit <- lmFit(v, design)
# contr <- makeContrasts(contrasts = paste0(levels(meta$Group)[2],"-",levels(meta$Group)[1]),
#                        levels=design)
# fit2 <- contrasts.fit(fit, contr)
# fit2 <- eBayes(fit2)

# res <- topTable(fit2, number = nrow(v$E), sort.by="P")
# # 저장
# write.csv(res, "pathway_limma_results.csv", row.names = TRUE)

# # volcano plot (pathway)
# res$logFC <- res$logFC
# res$negLogP <- -log10(res$P.Value)
# res$pathway <- rownames(res)

# library(ggplot2)
# ggplot(res, aes(x=logFC, y=negLogP)) +
#   geom_point(alpha=0.6) +
#   geom_text(data=head(res[order(res$P.Value),], 10), aes(label=pathway), size=3, vjust=1) +
#   theme_minimal() + labs(title="Pathway differential (limma)", x="log2FC", y="-log10(p)")

# # barplot top significant pathways
# topk <- head(res[order(res$P.Value),], 20)
# topk$pathway <- factor(rownames(topk), levels = rownames(topk))
# ggplot(topk, aes(x=pathway, y=logFC)) +
#   geom_bar(stat="identity") + coord_flip() +
#   theme_minimal() + labs(title="Top 20 differential pathways (by p-value)")
