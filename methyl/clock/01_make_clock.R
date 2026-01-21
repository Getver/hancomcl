library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(reshape2)
library(data.table)
library('ggsignif')
library(reshape)
library(stringr)
library(grid)
library(stringr)
library(limma)
library(tidyverse)
library(data.table)
library(ExperimentHubData)
library(ExperimentHub)
library(impute)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("glmnet")
# install.packages('glmnet')

source('/home/data/NTpair_1400/02_function.R')

library(glmnet)


get_experimenthub()
no
no
no
no
no
no
no

## after boruta pipeline(total 190,150 marker)
load('/home/data/data/beta_2024.rda')  # norm_by_bmiq  [1] 936990   2106
ls()

beta23 <- data.frame(fread('/home/data/data/beta_2023.txt', sep='\t', header=TRUE))
beta24 <- data.frame(norm_by_bmiq)

df <- cbind(beta23, beta24)
df$cpgs <- gsub('_.*$', '', row.names(df))
df_unique <- df %>% distinct(cpgs, .keep_all = TRUE)

sample <- read.csv('/home/data/data/sample_sheet.csv', header=TRUE, skip=6)

all(sample$Sample_Name == names(df_unique))  # 전체 순서 동일한지 확인
target <- intersect(sample$Sample_Name, names(df_unique))


# 936990 -> 930596  3826

inputBeta <- make_inputBeta(df_unique, df_unique$cpgs, target) # 2473 83060
all(row.names(inputBeta) == target) # 한번 더 확인
all(target == sample$Sample_Name)
inputAge <- sample$Age



###################################
## clock 생성
key = '01_935clock'
coef <- make_clock(inputBeta, inputAge, key)

### AGE 계산
coef_SM <- paste0('/home/data/candyCLOCK/age_glm_', key,'.txt')
coef_input <- read.table(coef_SM, sep='\t', header=TRUE, col.names=c('CpGmarker', 'CoefficientTraining'))


out <- go_analysis(df_unique, coef_input)

write.table(out, sep='\t', paste0('/home/data/candyCLOCK/DNAmAge_predict_', key, '.txt'), col.names=TRUE, row.names=FALSE, quote=FALSE)


# intersect(sumin$cpgs, coefHorvath$CpGmarker)

