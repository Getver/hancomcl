# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k
# https://partrita.github.io/posts/seurat-scRNAseq/         # pipeline
# https://m.blog.naver.com/ruins0408/222971721429           # annotation
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

library(Seurat)
# library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
package_version(R.version)
packageVersion("Seurat")

# SAMPLE <- args[1];
SAMPLE <- 'sc5p_v2_hs_PBMC_10k_5gex'
# SAMPLE <- 'pbmc8k'


dir.create(dirname(paste0("./", SAMPLE, "/plot/test.txt")), recursive = TRUE, showWarnings = FALSE)

data <- Read10X(data.dir = paste0("/home/data/",  SAMPLE, "/outs/filtered_feature_bc_matrix"))
raw_df <- CreateSeuratObject(counts = data, project = SAMPLE, min.cells = 3, min.features = 200)
raw_df[["percent.mt"]] <- PercentageFeatureSet(object = raw_df, pattern = "^MT-")

png(paste0("./", SAMPLE, "/plot/plot01_before_process.png"), width=1500, height=750)
VlnPlot(object = raw_df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


filtered_df <- subset(x = raw_df, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 50)

png(paste0("./", SAMPLE, "/plot/plot02_after_process.png"), width=1500, height=750)
VlnPlot(object = filtered_df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


normalized_df <- NormalizeData(object = filtered_df, normalization.method = "LogNormalize", scale.factor = 10000)
featured_df <- FindVariableFeatures(object = normalized_df, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = featured_df), 10)

png(paste0("./", SAMPLE, "/plot/plot03_highly_variable.png"), width=1500, height=750)
plot1 <- VariableFeaturePlot(object = featured_df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


all.genes <- rownames(x = featured_df)
scaled_df <- ScaleData(object = featured_df, features = all.genes)
PCA_df <- RunPCA(object = scaled_df, features = VariableFeatures(object = scaled_df))
# print(PCA_df[["pca"]], dims = 1:5, nfeatures = 5)


png(paste0("./", SAMPLE, "/plot/plot04_elbowplot.png"), width=600, height=600)
ElbowPlot(object = PCA_df)
dev.off()


neighbor_df <- FindNeighbors(object = PCA_df, dims = 1:10)
head(neighbor_df@graphs)
cluster_df <- FindClusters(object = neighbor_df, resolution = 0.5)
head(Idents(cluster_df), 5)
umap_df <- RunUMAP(object = cluster_df, dims = 1:10)

png(paste0("./", SAMPLE, "/plot/plot05_umap.png"), width=600, height=600)
DimPlot(object = umap_df, reduction = "umap")
dev.off()


dir.create(dirname(paste0("./", SAMPLE, "/cluster/test.txt")), recursive = TRUE, showWarnings = FALSE)
cluster0.markers <- FindMarkers(object = umap_df, ident.1 = 0, min.pct = 0.25)
write.table(cluster0.markers,file=paste0("./", SAMPLE, '/cluster/cluster0_marker.txt'),sep='\t',quote=FALSE)
cluster1.markers <- FindMarkers(object = umap_df, ident.1 = 1, min.pct = 0.25)
write.table(cluster1.markers,file=paste0("./", SAMPLE, '/cluster/cluster1_marker.txt'),sep='\t',quote=FALSE)
cluster2.markers <- FindMarkers(object = umap_df, ident.1 = 2, min.pct = 0.25)
write.table(cluster2.markers,file=paste0("./", SAMPLE, '/cluster/cluster2_marker.txt'),sep='\t',quote=FALSE)
cluster3.markers <- FindMarkers(object = umap_df, ident.1 = 3, min.pct = 0.25)
write.table(cluster3.markers,file=paste0("./", SAMPLE, '/cluster/cluster3_marker.txt'),sep='\t',quote=FALSE)
cluster4.markers <- FindMarkers(object = umap_df, ident.1 = 4, min.pct = 0.25)
write.table(cluster4.markers,file=paste0("./", SAMPLE, '/cluster/cluster4_marker.txt'),sep='\t',quote=FALSE)
cluster5.markers <- FindMarkers(object = umap_df, ident.1 = 5, min.pct = 0.25)
write.table(cluster5.markers,file=paste0("./", SAMPLE, '/cluster/cluster5_marker.txt'),sep='\t',quote=FALSE)

markers <- FindAllMarkers(object = umap_df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file=paste0("./", SAMPLE, '/cluster/markers.txt'),sep='\t',quote=FALSE)

tg <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
tgs <- tg$gene


png(paste0("./", SAMPLE, "/plot/plot06_vln.png"), width=1200, height=1200)
VlnPlot(object = umap_df, features = tgs)
dev.off()

png(paste0("./", SAMPLE, "/plot/plot07_umap.png"), width=1200, height=1200)
FeaturePlot(object = umap_df, features = tgs)
dev.off()

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
png(paste0("./", SAMPLE, "/plot/plot08_heatmap.png"), width=1200, height=1200)
DoHeatmap(object = umap_df, features = top10$gene) + NoLegend()
dev.off()


####################################

# https://jackbibby1.github.io/SCPA/






################################# 아래는 cell type annotation 방법

library(future)
library(SeuratDisk)
options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1GB

SCT_df <- SCTransform(filtered_df, verbose = FALSE) # 정규화
SCT_df <- RunPCA(SCT_df, features = VariableFeatures(object = SCT_df))
SCT_df <- RunUMAP(SCT_df, dims = 1:50)


reference_data <- LoadH5Seurat('/usr/local/src/reference/multi.h5seurat')
# reference_data <- read.csv('/usr/local/src/reference/Differential_expression_by_celltype.csv', header=T, sep=',')
# reference_df <- reference_data  # 예시로 동일 데이터 사용
# reference_df <- RunPCA(reference_df, features = VariableFeatures(object = reference_df)) # PCA 수행
# reference_df <- RunUMAP(reference_df, dims = 1:10)  # Reference 데이터에도 UMAP 적용


png(paste0("./", SAMPLE, "/plot/plot10_reference.png"), width=600, height=600)
DimPlot(object = reference_data, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
dev.off()



anchors <- FindTransferAnchors(
  reference = reference_data,
  query = SCT_df,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# head(reference_df@meta.data)

mapped_df <- MapQuery(
  anchorset = anchors,
  query = SCT_df,
  reference = reference_data,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)



png(paste0("./", SAMPLE, "/plot/plot09_predictedmap.png"), width=1500, height=750)
p1 = DimPlot(mapped_df, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(mapped_df, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2
dev.off()

# FeaturePlot(pbmc3k, features = c("pDC", "CD16 Mono", "Treg"),  reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))



# marker_list <- list()
# gene_number <- 10 # 10개의 유전자만 찾아낼때

# # for loop으로 리스트에 값 추가
# for (i in unique(Idents(umap_df))) {
#   marker_list[[paste0("cluster",i)]] <- markers %>% group_by(cluster) %>% top_n(gene_number, avg_log2FC) %>%
#       ungroup() %>% arrange(cluster, desc(avg_log2FC)) %>% filter(cluster == i) %>% .$gene
# }

# # 결과 출력(chatGPT에게 물어보기)
# marker_list



# 참조를 위해 이전 ID 클래스(클러스터 레이블)를 저장합니다.
# umap_df[["old.ident"]] <- Idents(object = umap_df)

# # 레이블 변경하기
# umap_df <- RenameIdents(
#     object = umap_df,
#     `0` = "neutrophils",
#     `1` = "CD4+ T cells",
#     `2` = "naive B cells",
#     `3` = "Plasma cells",
#     `4` = "monocytes",
#     `5` = "CD56bright NK cells",
#     `6` = "memory  B cells",
#     `7` = "CD8 T cells",
#     `8` = "CD56dim NK cells",
#     `9` = "platelet"
#     )


# p <- DimPlot(umap_df, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# options(repr.plot.width = 7, repr.plot.height = 7)

# table(Idents(umap_df))

