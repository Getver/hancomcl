library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(doParallel)

basedir <- '/home/DATA'

registerDoParallel(cores = 48)

library(ChAMP)
library(ChAMPdata)

targets <- read.metharray.sheet(basedir, pattern = "csv$", ignore.case = TRUE,
  recursive = TRUE, verbose = TRUE)

rgset <- read.metharray.exp(targets = targets, force = TRUE, extended = TRUE)
sampleNames(rgset) <- targets$Sample_Name

det_p <- detectionP(rgset, type = "m+u")

det_p_with_name <- cbind(rownames(det_p), det_p)
colnames(det_p_with_name)[1] <- "cpgs"

write.table(det_p_with_name, file = "/home/output/Minfi_detection-P_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

raw_mset <- preprocessRaw(rgset)
raw_beta <- getBeta(raw_mset)

raw_beta_with_name <- cbind(rownames(raw_beta), raw_beta)
colnames(raw_beta_with_name)[1] <- "cpgs"

write.table(raw_beta_with_name, file = "/home/output/raw_beta_value_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

norm_by_noob <- preprocessNoob(rgset, offset = 15, dyeCorr = TRUE,
  verbose = TRUE, dyeMethod = c("single"))

norm_by_bmiq <- champ.norm(beta = getBeta(norm_by_noob), mset = norm_by_noob,
  method = "BMIQ", cores = 1, plotBMIQ = FALSE, arraytype='EPICv1')

norm_by_bmiq_with_name <- cbind(rownames(norm_by_bmiq), norm_by_bmiq)
colnames(norm_by_bmiq_with_name)[1] <- "cpgs"

write.table(norm_by_bmiq_with_name, file = "/home/output/two-way_normalized_beta_value_all_CpGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
