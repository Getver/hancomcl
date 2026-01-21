# library(minfi)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(reshape2)
library(data.table)
library(pheatmap)
library('ggsignif')
library(reshape)
library(stringr)
library(grid)
library(stringr)
library(limma)
library(glmnet)
library(tidyverse)


sample <- read.table('/disk0/sm/methyl/03_pollution/total/work/sample_sheet_re.txt', sep='\t', header=T)


## TEST1_1056

set.seed(1234)

# 각 나이별 22명씩 추출
age <- data.frame(table(sample$Age))[1:48,]
sample <- sample %>% filter(Age < 88)
split_df <- split(sample, sample$Age)
filtered_list <- lapply(split_df, function(x) x[sample(nrow(x), 22), ])
result <- do.call(rbind, filtered_list)

# 2023년 데이터 호출
load('/disk0/sm/methyl/03_pollution/total/work/B_beta.Rdata')   # norm_by_bmiq
target <- result[(result$Pool_ID == '2023'),]
b2023 <- norm_by_bmiq[, target$Sample_Name]

# 2024년 데이터 호출
load('/disk0/sm/methyl/03_pollution/total/B_beta_nor.rda')   # norm_by_bmiq
target <- result[(result$Pool_ID == '2024'),]
b2024 <- norm_by_bmiq[, target$Sample_Name]

beta <- cbind(b2023, b2024)

row.names(beta) <- sub("_.+$", "", row.names(beta))
e4 <- fread('/disk0/sm/methyl/TEST/1010_regression_3/cpgs_450K.txt', header=TRUE, sep = '\t')$Illumina
epic1 <- fread('/disk0/sm/methyl/TEST/1010_regression_3/EPIC1_loci.txt', header=TRUE, sep = '\t')$cpgs
epic2 <- sub("_.+$", "", fread('/disk0/sm/methyl/TEST/1010_regression_3/cpgs_epic2.txt', header=TRUE, sep = '\t')$cpgs)
e <- intersect(intersect(epic1, epic2), e4)
ebeta <- beta[e, ]
tbeta <- t(ebeta)    # [, 1:20000]
tmpbeta <- data.frame(tbeta)

result$Age <- as.numeric(result$Age)
pheno <- (result$Age - 20)/(21)

key='TEST1_1056'

glmnet.Training.CV <- cv.glmnet(as.matrix(tmpbeta), pheno, nfolds=10, alpha=0.5,family="gaussian")
lambda.glmnet.Training <- glmnet.Training.CV$lambda.min
glmnet.Training <- glmnet(as.matrix(tmpbeta), pheno, family="gaussian", alpha=0.5, nlambda=10, lambda=lambda.glmnet.Training)# 0.042, lambda=lambda.glmnet.Training

co <- coef(glmnet.Training, s=lambda.glmnet.Training)
cpg <- c()
value <- c()
for (i in 1:length(co)){
    # print(i)
    if (co[i] != 0){
        cpg <- c(cpg, rownames(co)[i])
        value <- c(value, co[i])
    }
}

table <- data.frame(cbind(cpg, value))
write.table(table, paste0('/disk0/sm/methyl/03_pollution/total/clock/age_glm_', key,'.txt'), sep='\t', col.names=T, row.names=F, quote=F)





## TEST2_noob
load('/disk0/sm/methyl/03_pollution/total/B_beta_noob.rda') # norm_by_noob


beta <- cbind(norm_by_noob)

row.names(beta) <- sub("_.+$", "", row.names(beta))
ebeta <- beta[e, ]
tbeta <- t(ebeta)    # [, 1:20000]
tmpbeta <- data.frame(tbeta)

target <- result[(result$Pool_ID == '2024'),]
target$Age <- as.numeric(target$Age)
pheno <- (target$Age - 20)/(21)


key='TEST2_noob'


glmnet.Training.CV <- cv.glmnet(as.matrix(tmpbeta), pheno, nfolds=10, alpha=0.5,family="gaussian")
lambda.glmnet.Training <- glmnet.Training.CV$lambda.min
glmnet.Training <- glmnet(as.matrix(tmpbeta), pheno, family="gaussian", alpha=0.5, nlambda=10, lambda=lambda.glmnet.Training)# 0.042, lambda=lambda.glmnet.Training

co <- coef(glmnet.Training, s=lambda.glmnet.Training)
cpg <- c()
value <- c()
for (i in 1:length(co)){
    # print(i)
    if (co[i] != 0){
        cpg <- c(cpg, rownames(co)[i])
        value <- c(value, co[i])
    }
}

table <- data.frame(cbind(cpg, value))
write.table(table, paste0('/disk0/sm/methyl/03_pollution/total/clock/age_glm_', key,'.txt'), sep='\t', col.names=T, row.names=F, quote=F)






## bedtools

load('/disk0/sm/methyl/03_pollution/total/clock/B_beta_37K.Rdata')  # ebeta


meth_annotation <- function(path) {
  EPIC2 <- read.csv(path, sep=',', header=FALSE)
  names(EPIC2) <- EPIC2[8,]
  e <- EPIC2[9:nrow(EPIC2),]
  e2 <- e[c('IlmnID', 'CHR', 'MAPINFO', 'UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'UCSC_RefGene_Name')]
  e2$UCSC_RefGene_Group <- sub(";.+$", "", e2$UCSC_RefGene_Group)
  e2$UCSC_RefGene_Group <- sub("_.+$", "", e2$UCSC_RefGene_Group)
  e2$UCSC_RefGene_Name <- sub(";.+$", "", e2$UCSC_RefGene_Name)
  e2$Relation_to_UCSC_CpG_Island <- sub(";.+$", "", e2$Relation_to_UCSC_CpG_Island)
  return(e2)
}

# K450 <- meth_annotation('/disk0/sm/methyl/mchip/HM450.csv')
# epic1 <- meth_annotation('/disk0/sm/methyl/mchip/EPICv1.csv')
epic2 <- meth_annotation('/disk0/sm/methyl/mchip/EPICv2.csv')
epic2$end <- as.numeric(epic2$MAPINFO) + 1
bed <- epic2[,c('CHR', 'MAPINFO', 'end', 'IlmnID')]
bed$IlmnID <- sub("_.+$", "", bed$IlmnID)


bed <- bed %>% filter(IlmnID %in% row.names(data.frame(ebeta)))
sed <- bed[order(bed$CHR, bed$MAPINFO), ]
sed <- sed %>% filter(CHR != 'chr0')

write.table(sed, '/disk0/sm/37K.bed', sep='\t', row.names=F, col.names=F, quote=F)

system('bedtools intersect -a /disk0/sm/37K.bed -b /disk0/sm/methyl/03_pollution/total/clock/CpGisland.txt -wa > /disk0/sm/methyl/03_pollution/total/clock/CpGisland_37K.txt')
island_cpg <- read.table('/disk0/sm/methyl/03_pollution/total/clock/CpGisland_37K.txt', sep='\t', header=F) # 116175


# promoter <- read.table('/disk0/kb/temp/Gencode_v47_promoter.txt', header=T, sep='\t')
system('grep -v "_" /disk0/kb/temp/Gencode_v47_promoter.txt | cut -f -3 | grep -v "-" > /disk0/sm/methyl/03_pollution/total/clock/promoter_gencode_V47.txt')
system('bedtools intersect -a /disk0/sm/37K.bed -b /disk0/sm/methyl/03_pollution/total/clock/promoter_gencode_V47.txt -wa > /disk0/sm/methyl/03_pollution/total/clock/promoter_37K.txt')
promoter_cpg <- read.table('/disk0/sm/methyl/03_pollution/total/clock/promoter_37K.txt', sep='\t', header=F) # 116175

e2 <- epic2 %>% filter(UCSC_RefGene_Group == 'exon')
e2$IlmnID <- sub("_.+$", "", e2$IlmnID)

y <- unique(c(promoter_cpg$V4, island_cpg$V4, e2$IlmnID))


## 
dmp <- read.table('./DMP_Minfi_dmpFinder.txt', sep='\t', header=T)
qmp <- dmp %>% filter((qval < 0.00000001))
mp <-  dmp %>% filter((abs(t) > 5))


y <- unique(sub("_.+$", "", c(rownames(mp), rownames(qmp)))) 





