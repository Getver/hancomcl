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

# Specialized Methylation Aging Determiner for Koreans

# load("B_beta_nor.rda")
# e9 <- norm_by_bmiq
# row.names(e9) <- sub("_.+$", "", row.names(norm_by_bmiq))
# e4 <- fread('/disk0/sm/methyl/TEST/1010_regression_3/cpgs_450K.txt', header=TRUE, sep = '\t')$Illumina
# epic1 <- fread('/disk0/sm/methyl/TEST/1010_regression_3/EPIC1_loci.txt', header=TRUE, sep = '\t')$cpgs
# epic2 <- sub("_.+$", "", fread('/disk0/sm/methyl/TEST/1010_regression_3/cpgs_epic2.txt', header=TRUE, sep = '\t')$cpgs)
# e <- intersect(intersect(epic1, epic2), e4)
# ebeta <- e9[e, ]
# save(ebeta, file = "./clock/B_beta_37K.Rdata")

load("./clock/B_beta_37K.Rdata")


# elastic_glm <- function(ebeta, number, key, group, no) {



# sample <- read.table('./work/sample_sheet_new.csv', sep=",", fill = TRUE, col.names=c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position', 'Group', 'Age', 'Gender'))[8:2011, ]
sample <- read.table('/disk0/sm/methyl/03_pollution/total/sample_sheet.csv', sep=",", fill = TRUE, col.names=c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position', 'Group', 'Gender', 'Age'))[8:2113, ]
target <- sample[(sample$Pool_ID == '2024'),]
target$Age <- as.numeric(target$Age)
pheno <- (target$Age - 20)/(21)
beta <- ebeta[, target$Sample_Name]
# beta <- norm_by_bmiq[, target$Sample_Name]
tbeta <- t(beta)    # [, 1:20000]

key='dmp_qval_55'

set.seed(1234)
# set.seed(1234)
tmpbeta <- data.frame(tbeta)

# no <- read.table('./clock/no.txt', sep='\t')
# for (i in no){
#     tmpbeta[,i] <- NULL
# } 

# yes <- read.table('./clock/yes.txt', sep='\t')
# tmpbeta <- tmpbeta[,intersect(y, names(tmpbeta))]
# y <- intersect(etarget$IlmnID, names(tmpbeta))    # 56462

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

# }
# save.image(file = "./clock/S2106_elasticnet.RData")

# e3 <- e2 %>% filter(Relation_to_UCSC_CpG_Island != '') %>% filter(UCSC_RefGene_Group != '')
# e3$id <- sub("_.+$", "", e3$IlmnID)
# yes <- intersect(e3$id, names(tmpbeta))
# tmpbeta <- tmpbeta[,yes]


# e4
# epic1

# a$cpg <- sub("_.+$", "", a$cpg)
# a1 <- intersect(a$cpg, e4)
# a2 <- intersect(a$cpg, epic1)

# b <- a %>% filter(!(cpg %in% a1)) %>% filter(!(cpg %in% a2))


## Recursive Feature Elimination(RFE)

# library(caret)

# glmnet.Training.CV <- cv.glmnet(as.matrix(tmpbeta), pheno, nfolds=10, alpha=0.5,family="gaussian")
# lambda.glmnet.Training <- glmnet.Training.CV$lambda.min

# control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

# rfe_results  <- rfe(as.matrix(tmpbeta), pheno, sizes = c(1:ncol(tmpbeta)), rfeControl = control, method='glmnet', tuneGrid = expand.grid(alpha = 0.5, lambda=lambda.glmnet.Training))




## PCA


# fit <- cv.glmnet(as.matrix(tmpbeta), pheno, nfolds=10, alpha=0.5,family="gaussian")
# best_lambda <- fit$lambda.min

# coef(fit, s = "lambda.min")

# selected_vars <- data.frame(as.matrix(coef(fit, s=best_lambda)))
# non_zero_vars <- selected_vars %>% filter(s1 != 0)

# X_selected <- tmpbeta %>% select(row.names(non_zero_vars)[2:nrow(non_zero_vars)])
# pca_result <- prcomp(X_selected, center = TRUE, scale. = TRUE)

# summary(pca_result)

# X_pca <- pca_result$x
# fit_pca <- cv.glmnet(as.matrix(X_pca), pheno, nfolds=10, alpha=0.5,family="gaussian")

# coco <- coef(fit_pca, s = "lambda.min")
# noco <- data.frame(as.matrix(coco)) %>% filter(s1 != 0)

# # head(pca_result$rotation)

# result <- pca_result$rotation
# write.table(result, './rotation.txt', sep='\t', quote=F)


# png(filename="myplot.png",width=300,height=600,unit="px",bg="transparent")
# plot(pca_result, type="l")
# dev.off()






## horvath랑 sm annotation 비교하기(for go)

horvath <- read.csv('/disk0/sm/methyl/TEST/1031/Horvath_gene.txt', sep='\t', header=T)
hcpg <- horvath$CpGmarker

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

K450 <- meth_annotation('/disk0/sm/methyl/mchip/HM450.csv')
horvath_df <- K450 %>% filter(IlmnID %in% hcpg)
horvath_df$group <- 'Horvath'
epic1 <- meth_annotation('/disk0/sm/methyl/mchip/EPICv1.csv')
epic2 <- meth_annotation('/disk0/sm/methyl/mchip/EPICv2.csv')
epic2$IlmnID <- sub("_.+$", "", epic2$IlmnID)
sm <- read.csv('/disk0/sm/methyl/TEST/1031/train/train3_age_glm.txt', sep='\t', header=T)

sm_df <- epic2 %>% filter(IlmnID %in% sm$cpg)
sm_df$group <- 'SM'

# both <- intersect(horvath_df$IlmnID, sm_df$IlmnID)\
# horvath_df %>% filter(IlmnID %in% both)

horvath_t <- data.frame(table(horvath_df$UCSC_RefGene_Name))
sm_t <- data.frame(table(sm_df$UCSC_RefGene_Name))

horvath_t$group <- 'Horvath'
sm_t$group <- 'SM'

gene_count <- rbind(horvath_t, sm_t)

# intersect(horvath_df$UCSC_RefGene_Name, sm_df$UCSC_RefGene_Name)    # 341 377   16

write.table(gene_count, '/disk0/sm/methyl/TEST/1031/data/gene_count.txt', sep='\t', row.names=F, col.names=T, quote=F)


nrow(epic2 %>% filter(UCSC_RefGene_Name == 'ELOVL2'))   #ELOVL2 FHL2
nrow(epic2 %>% filter(UCSC_RefGene_Name == 'ELOVL2-AS1'))


## G456 G678 겹치는 마커 비교

G45 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G1_45.txt', sep='\t', header=T)
G46 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G1_46.txt', sep='\t', header=T)
G56 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G1_56.txt', sep='\t', header=T)
G67 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G11_67.txt', sep='\t', header=T)
G68 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G11_68.txt', sep='\t', header=T)
G78 <- read.csv('/disk0/sm/methyl/TEST/1031/Age/age_glm_G11_78.txt', sep='\t', header=T)


# FG45 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG1_45.txt', sep='\t', header=T)
# FG46 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG1_46.txt', sep='\t', header=T)
# FG56 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG1_56.txt', sep='\t', header=T)
# FG67 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG11_67.txt', sep='\t', header=T)
# FG68 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG11_68.txt', sep='\t', header=T)
# FG78 <- read.table('/disk0/sm/methyl/TEST/1017_age/pvalFilter/001/DMP_Minfi_dmpFinder_FG11_78.txt', sep='\t', header=T)

# count_cpg <- function(data, Fdata){
#     Fdata$V1 <- sub("_.+$", "", Fdata$V1)
#     length(intersect(data$cpg, Fdata$V1))
# }

# count_cpg(G45, FG45)
# count_cpg(G46, FG46)
# count_cpg(G56, FG56)
# count_cpg(G67, FG67)
# count_cpg(G68, FG68)
# count_cpg(G78, FG78)

# lapply(c(G45, G46, G56, G67, G68, G78), intersect())
# all <- intersect(intersect(intersect(intersect(intersect(G45$cpg, G46$cpg), G56$cpg), G67$cpg), G68$cpg), G78$cpg)  # 6

G <- list(G45, G46, G56, G67, G68, G78)
GG <- Reduce(function(x, y) merge(x, y, by='cpg', all=T), G)
names(GG) <- c('cpg', 'G45', 'G46', 'G56', 'G67', 'G68', 'G78')

write.table(GG, '/disk0/sm/methyl/TEST/1031/data/G45678_cpg.txt', sep='\t', row.names=F, col.names=T, quote=F)

al <- list(G45, G46, G56, G67, G68, G78)
all <- al %>% reduce(inner_join, by='cpg')


length(intersect(G67$cpg, G78$cpg))



# data = reduce(lambda x,y: pd.merge(x,y, on=['TargetID'], how='inner'), [data1, data2, data3, data4, data5])




## 3만 4천개 뽑기, chr보기

G40 <- read.csv('/disk0/sm/methyl/TEST/1023_gender/DMP_Minfi_dmpFinder_FG40.txt', sep='\t', header=T)
G50 <- read.csv('/disk0/sm/methyl/TEST/1023_gender/DMP_Minfi_dmpFinder_FG50.txt', sep='\t', header=T)
G60 <- read.csv('/disk0/sm/methyl/TEST/1023_gender/DMP_Minfi_dmpFinder_FG60.txt', sep='\t', header=T)

FG40 <- epic2 %>% filter(IlmnID %in% G40$V1)
FG50 <- epic2 %>% filter(IlmnID %in% G50$V1)
FG60 <- epic2 %>% filter(IlmnID %in% G60$V1)

al <- list(FG40, FG50, FG60)
all <- al %>% reduce(inner_join, by='IlmnID')


only <- all %>% filter((CHR != 'chrX') & (CHR != 'chrY'))



gbeta <- norm_by_bmiq[all$IlmnID,]


sample <- read.table('/disk0/sm/methyl/TEST/1031/work/sample_sheet_new.csv', sep=",", fill = TRUE, col.names=c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position', 'Group', 'Age', 'Gender'))[8:3731, ]
targets <- sample[order(sample$Sample_Name), ]
target <- targets %>% filter(Pool_ID == 2024)
target$pheno <- 0
target[target$Gender == 'M', 'pheno'] = 1

pheno <- target$pheno
beta <- gbeta[, target$Sample_Name]
tbeta <- t(beta)    # [, 1:20000]
glmnet.Training.CV <- cv.glmnet(tbeta, pheno, nfolds=10, alpha=0.5,family="gaussian")
lambda.glmnet.Training <- glmnet.Training.CV$lambda.min
glmnet.Training <- glmnet(tbeta, pheno, family="gaussian", alpha=0.5, nlambda=100)
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

write.table(table, paste0('/disk0/sm/methyl/TEST/1031/Age/GENDER_elasticnet.txt'), sep='\t', col.names=T, row.names=F, quote=F)



intersect(all$IlmnID, table$cpg)

