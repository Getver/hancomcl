library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(RColorBrewer)

print(warnings())

arguments <- commandArgs(trailingOnly = TRUE)
cwd <- arguments[1]
platform <- arguments[2]

setwd(cwd)

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

basedir <- cwd

targets <- read.metharray.sheet(basedir, pattern = "csv$", ignore.case = TRUE,
  recursive = TRUE, verbose = TRUE)

rgset <- read.metharray.exp(targets = targets, force = TRUE, extended = TRUE)

sampleNames(rgset) <- targets$Sample_Name

if (platform == "EPICv2") {
  rgset@annotation <- c(array = "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38")
}

loci <- getAnnotation(rgset)
df_loci <- data.frame(loci)

if (platform == "EPICv2") {
  loci <- df_loci[c("chr", "pos", "strand", "Name", "Probe_rs", "CpG_rs",
    "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "Methyl450_Loci",
    "Methyl27_Loci", "EPICv1_Loci")]
} else {
  loci <- df_loci[c("chr", "pos", "strand", "Name", "Probe_rs", "CpG_rs",
    "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "Methyl450_Loci",
    "Methyl27_Loci")]
}

write.table(loci, file = "loci.txt", sep = "\t", col.names = TRUE,
  row.names = FALSE, quote = FALSE)

det_p <- detectionP(rgset, type = "m+u")

det_p_with_name <- cbind(rownames(det_p), det_p)
colnames(det_p_with_name)[1] <- "cpgs"

write.table(det_p_with_name, file = "Minfi_detection-P_all_CpGs.txt",
  sep = "\t", row.names = FALSE, quote = FALSE)

raw_mset <- preprocessRaw(rgset)
raw_beta <- getBeta(raw_mset)

raw_beta_with_name <- cbind(rownames(raw_beta), raw_beta)
colnames(raw_beta_with_name)[1] <- "cpgs"

write.table(raw_beta_with_name, file = "raw_beta_value_all_CpGs.txt",
  sep = "\t", row.names = FALSE, quote = FALSE)

norm_by_noob <- preprocessNoob(rgset, offset = 15, dyeCorr = TRUE,
  verbose = TRUE, dyeMethod = c("single"))

# normalize by BMIQ in ChAMP followed by Noob nomalization
norm_by_bmiq <- champ.norm(beta = getBeta(norm_by_noob), mset = norm_by_noob,
  method = "BMIQ", cores = 48, arraytype = platform, plotBMIQ = FALSE)

norm_by_bmiq_with_name <- cbind(rownames(norm_by_bmiq), norm_by_bmiq)
colnames(norm_by_bmiq_with_name)[1] <- "cpgs"

write.table(norm_by_bmiq_with_name,
  file = "two-way_normalized_beta_value_all_CpGs.txt",
  sep = "\t", row.names = FALSE, quote = FALSE)


# retrieve cpgs for filter by detection p value
p_failed <- det_p > 0.05

p_removed <- names(which(rowMeans(p_failed) > 0, TRUE))

remove_probes <- match(p_removed, rownames(norm_by_noob)) %>% unique %>% na.omit
filtered_norm_by_noob <- norm_by_noob[-remove_probes, ] # MethylSet