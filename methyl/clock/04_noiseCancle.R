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


load("./work/B_beta.Rdata")
load("./work/B_beta_37K.Rdata")



sample <- read.table('./work/sample_sheet_new.csv', sep=",", fill = TRUE, col.names=c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position', 'Group', 'Age', 'Gender'))[8:3731, ]
target <- sample[(sample$Pool_ID == '2024'),]
target$Age <- as.numeric(target$Age)
pheno <- (target$Age - 20)/(21)
beta <- ebeta[, target$Sample_Name]
tbeta <- t(beta)    # [, 1:20000]

key='S2106_n015_2000'

set.seed(1234)
tmpbeta <- data.frame(tbeta)


correl_total <- read.table('/disk0/kb/temp/meth/2023_meth_array/correlation_1720.txt', sep='\t', header=T)  # correlation 분석 결과(cpg-age)
correl_pval <- correl_total %>% filter(p_adjusted < 0.05)

corr <- correl_pval %>% filter(abs(cor_results) >= 0.15)   # 122024
el <- correl_pval[order(abs(correl_pval$cor_results)),][1:2000,]
correl <- rbind(corr, el)
correl$cpgs <- sub("_.+$", "", correl$cpgs)
cpg <- intersect(correl$cpgs, names(tmpbeta))  # 75118
tmpbeta <- tmpbeta[,cpg]

# yes <- read.table('/disk0/sm/methyl/TEST/1031_markerDel/yes.txt', sep='\t')
# tmpbeta <- tmpbeta[,yes$V1]

# no <- read.table('./no.txt', sep='\t')
# for (i in no){
#     tmpbeta[,i] <- NULL
# }



glmnet.Training.CV <- cv.glmnet(as.matrix(tmpbeta), pheno, nfolds=10, alpha=0.5,family="gaussian")
lambda.glmnet.Training <- glmnet.Training.CV$lambda.min
glmnet.Training <- glmnet(as.matrix(tmpbeta), pheno, family="gaussian", alpha=0.5, nlambda=100, lambda=lambda.glmnet.Training)# 

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
write.table(table, paste0('./', key,'.txt'), sep='\t', col.names=T, row.names=F, quote=F)

# }
