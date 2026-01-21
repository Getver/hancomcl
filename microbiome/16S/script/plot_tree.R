library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(data.table)

# 데이터 불러오기
a <- data.frame(fread('abundance/abundance.tsv', skip = 1))
tax <- data.frame(fread('taxonomy/taxonomy.tsv'))

# species 추출
tax$species <- sub(".*s__([^;]+).*", "\\1", tax$Taxon)

# strain 추출 (있다면)
tax$strain <- sub(".*t__([^;]+).*", "\\1", tax$Taxon)

# species + strain 합치기 (strain 없으면 species만)
tax$species_strain <- ifelse(tax$strain != "", paste0(tax$species, " (", tax$strain, ")"), tax$species)

# 매핑 (OTU ID → species_strain)
a$name <- tax$species_strain[match(a$X.OTU.ID, tax$Feature.ID)]

# OTU ID 고유하게 rownames 설정
rownames(a) <- make.unique(a$X.OTU.ID)

# 평균 abundance 계산
abundance_avg <- rowMeans(a[, sapply(a, is.numeric)])
mean_abundance <- mean(abundance_avg)

# 평균보다 높은 ASV 선택
top_asvs <- names(abundance_avg[abundance_avg > mean_abundance])

# tree 불러오기
tree <- read.tree("tree/tree.nwk")

# tree에서 top ASV만 선택
tree_sub <- drop.tip(tree, setdiff(tree$tip.label, top_asvs))

# OTU ID → species_strain 매핑
id_to_species <- setNames(tax$species_strain, tax$Feature.ID)

# tip.label을 species_strain으로 변경 (NA는 그대로 유지)
tree_sub$tip.label <- sapply(tree_sub$tip.label, function(x) {
  if (!is.na(id_to_species[x]) && id_to_species[x] != "") {
    id_to_species[x]
  } else {
    x
  }
})

# tree 그리기 (phylogram)
png("tree/tree_plot.png", width=2000, height=2000, res=300)
plot(tree_sub, type="phylogram", cex=0.6, tip.color="#848498", edge.width=1)
dev.off()

# tree 그리기 (circular)
png("tree/tree_plot_circular.png", width=2000, height=2000, res=300)
plot(tree_sub, type="fan", cex=0.6, tip.color="#848498", edge.width=1)
dev.off()
