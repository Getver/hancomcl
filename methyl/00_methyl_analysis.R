# docker run --rm -it -v ${PWD}:/home/data sm:meth_0.1

library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(doParallel)

registerDoParallel(cores=48)
print(warnings())

# platform <- "450K"
# platform <- "EPIC"
platform <- "EPICv2"

# setwd(cwd)

original_path <- .libPaths()

if (platform == "EPICv2") {
  ver2_path <- c(original_path[2], original_path[1])
  .libPaths(ver2_path)
  library(ChAMP)
  library(ChAMPdata)
  data(AnnoEPICv2)
} else {
  library(ChAMP)
  library(ChAMPdata)
}

basedir <- '/home/data'

targets <- read.metharray.sheet(basedir, pattern = "sample_sheet.csv$", ignore.case = TRUE,
  recursive = TRUE, verbose = TRUE)

# targets <- read.csv('/home/data/original_data/05_Sample_Sheet.csv', header=T)


## 샘플 빼기
# sample <- read.table('/home/data/data/pval_minfi.txt', sep='\t', header=T) %>% filter(check == 'FALSE')
# targets <- targets %>% filter(!Sample_Name %in% sample$X)



rgset <- read.metharray.exp(targets = targets, force = TRUE, extended = TRUE)
save(rgset, file = "B_rgset.rda")

# save.image(file = "sesame_rgset.RData")

sampleNames(rgset) <- targets$Sample_Name

if (platform == "EPICv2") {
  rgset@annotation <- c(array = "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38")
}

# loci <- getAnnotation(rgset)
# df_loci <- data.frame(loci)

# if (platform == "EPICv2") {
#   loci <- df_loci[c("chr", "pos", "strand", "Name", "Probe_rs", "CpG_rs",
#     "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "Methyl450_Loci",
#     "Methyl27_Loci", "EPICv1_Loci")]
# } else {
#   loci <- df_loci[c("chr", "pos", "strand", "Name", "Probe_rs", "CpG_rs",
#     "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "Methyl450_Loci",
#     "Methyl27_Loci")]
# }

# write.table(loci, file = "loci.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

det_p <- detectionP(rgset, type = "m+u")

det_p_with_name <- cbind(rownames(det_p), det_p)
colnames(det_p_with_name)[1] <- "cpgs"

write.table(det_p_with_name, file = "Minfi_detection-P_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

raw_mset <- preprocessRaw(rgset)
raw_beta <- getBeta(raw_mset)

raw_beta_with_name <- cbind(rownames(raw_beta), raw_beta)
colnames(raw_beta_with_name)[1] <- "cpgs"

write.table(raw_beta_with_name, file = "raw_beta_value_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# save(raw_beta_with_name, file = "B_beta_raw.rda")
# save.image(file = "raw_nor.RData")

norm_by_noob <- preprocessNoob(rgset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod = c("single"))
norm_by_noobbeta <- getBeta(norm_by_noob, type = "Illumina", offset = 0, betaThreshold = 0)


norm_by_noob_with_name <- cbind(rownames(norm_by_noobbeta), norm_by_noobbeta)
colnames(norm_by_noob_with_name)[1] <- "cpgs"

write.table(norm_by_noob_with_name, file = "Noob_beta_value_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)




# save(norm_by_noob, file = "B_beta_noob.rda")
# normalize by BMIQ in ChAMP followed by Noob nomalization # 이게 champ.DMR 로 들어감
norm_by_bmiq <- champ.norm(beta = getBeta(norm_by_noob), mset = norm_by_noob, method = "BMIQ", cores = 1, arraytype = platform, plotBMIQ = FALSE)
save(norm_by_bmiq, file = "B_beta_nor.rda")

norm_by_bmiq_with_name <- cbind(rownames(norm_by_bmiq), norm_by_bmiq)
colnames(norm_by_bmiq_with_name)[1] <- "cpgs"

write.table(norm_by_bmiq_with_name, file = "two-way_normalized_beta_value_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)



##################



Mset <- preprocessNoob(rgset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod=c("single"))

# save.image(file = "noob_nor_120.RData")
# load( "/home/data/noob_nor.RData" )


## 샘플 이름 바꾸기(KoBB id 변경)
# target <- data.frame(targets)
# target$ID <- gsub('-', '.', target$ID)
# data <- data.frame(pval)
# row.names(data) <- data$cpgs
# data$cpgs <- NULL

# col <- c()
# for (i in names(data)){
#   name <- target %>% filter(ID == i) %>% select(Sample_Name)
#   col <- c(col, name$Sample_Name)
# }
# names(data) <- col


########

Rset <- ratioConvert(Mset, what = "both", keepCN = TRUE)
GRset <- mapToGenome(Rset)
pheno <- as.character(targets$Group)
# pheno <- as.character(targets$Sample_Group)
designMatrix <- model.matrix(~pheno)

registerDoParallel(cores=4)

# subset_sumin <- subsetByLoci(rgset, includeLoci = rownames(raw_beta)[1:100], excludeLoci = NULL, keepControls = TRUE, keepSnps = TRUE)

# loci <- getAnnotation(rgset)
# df_loci <- data.frame(loci)

# minfi로 DMP 구하기
dmp_minfi <- dmpFinder(norm_by_bmiq, pheno, type = "continuous", qCutoff = 1)
write.table(dmp_minfi, file="DMP_Minfi_dmpFinder.txt", sep="\t", row.names=TRUE, quote=FALSE)
# dmp_minfi <- dmpFinder2(norm_by_bmiq, pheno, type = "categorical", qCutoff = 1)
# write.table(dmp_minfi, file="DMP_Minfi_dmpFinder2.txt", sep="\t", row.names=TRUE, quote=FALSE)

## minfi로 DMR 구하기
dmr_minfi <- bumphunter(GRset, design = designMatrix, cutoff = 0, B=0, type="Beta")
write.table(dmr_minfi$table, file="DMR_Minfi_bumphunter.txt", sep="\t", row.names=FALSE, quote=FALSE)


## champ로 DMP 구하기 ## 에러가 나네 ## 안나네
dmp_champ <- champ.DMP(norm_by_bmiq, pheno, arraytype='EPICv2', adjust.method='BH')
write.table(dmp_champ, file = "DMP_ChAMP.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## champ로 DMR 구하기
dmr_champ <- champ.DMR(norm_by_bmiq, pheno, arraytype='EPICv2', method='Bumphunter')
write.table(dmr_champ$BumphunterDMR, file = "DMR_ChAMP.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# DMR.GUI(DMR=dmr, beta=norm_by_bmiq, pheno=pheno, compare.group=NULL, arraytype='EPICv2')

# gender <- getSex(object = GRset, cutoff = -2)
# gender$gender <- pheno
# write.table(gender, file = "gender.txt", sep = "\t", quote = FALSE)

####################

# save.image(file = "noob_nor.RData")
# load( "/home/data/noob_nor.RData" )

# retrieve cpgs for filter by detection p value
# p_failed <- det_p > 0.05

# p_removed <- names(which(rowMeans(p_failed) > 0, TRUE))

# remove_probes <- match(p_removed, rownames(norm_by_noob)) %>% unique %>% na.omit
# filtered_norm_by_noob <- norm_by_noob[-remove_probes, ] # MethylSet





library(methylclockData)  # BiocManager::install("methylclock")
library(methylclock)
# library(Biobase)
# library(tibble)
# library(impute)
# library(ggplot2)
# library(ggpmisc)
# library(GEOquery)
library(data.table)

# MethylationData <- fread('./two-way_normalized_beta_value_all_CpGs.txt', header=TRUE, sep = '\t')
# MethylationData <- data.table(norm_by_bmiq_with_name)
MethylationData <- read.table('./two-way_normalized_beta_value_all_CpGs.txt', header=TRUE, sep = '\t')
MethylationData$cpgs <- sub("_.+$", "", MethylationData$cpgs)
# head(MethylationData$cpgs)
# a <- data.table(MethylationData)
# row.names(a) <- a$cpgs
# b <- a[2:length(a)]
checkClocks(MethylationData)

age.example55 <- DNAmAge(MethylationData)
# age.example55
write.table(age.example55, file='./DNAmAge_Result.txt', sep='\t', row.names = F, quote=F)

# age.example55.2 <- DNAmGA(MethylationData)
# # age.example55.2
# write.table(age.example55.2, file='./DNAmGA_Result.txt', sep='\t', row.names = F, quote=F)





# pvalue1 <- fread('/disk0/sm/methyl/03_pollution/total/2004/Minfi_detection-P_all_CpGs.txt', sep='\t', header=T)
# pvalue2 <- read.table('/disk0/sm/methyl/03_pollution/total/Minfi_detection-P_all_CpGs.txt', sep='\t', header=T)

# beta1 <- fread('/disk0/sm/methyl/03_pollution/total/2004/two-way_normalized_beta_value_all_CpGs.txt', sep='\t', header=T)
# beta2 <- read.table('/disk0/sm/methyl/03_pollution/total/two-way_normalized_beta_value_all_CpGs.txt', sep='\t', header=T)

# raw1 <- fread('/disk0/sm/methyl/03_pollution/total/2004/raw_beta_value_all_CpGs.txt', sep='\t', header=T)
# raw2 <- read.table('/disk0/sm/methyl/03_pollution/total/raw_beta_value_all_CpGs.txt', sep='\t', header=T)

# pvalue <- merge(pvalue1, pvalue2, by='cpgs')
# beta <- merge(beta1, beta2, by='cpgs')
# raw <- merge(raw1, raw2, by='cpgs')

# write.table(pvalue, '/disk0/sm/methyl/03_pollution/total/2100_pvalue.txt', sep='\t', quote=F, row.names=F)
# write.table(beta, '/disk0/sm/methyl/03_pollution/total/2100_beta.txt', sep='\t', quote=F, row.names=F)
# write.table(raw, '/disk0/sm/methyl/03_pollution/total/2100_raw.txt', sep='\t', quote=F, row.names=F)


