library("edgeR")
# args<-commandArgs(trailingOnly=TRUE);
# Name <- args[1];
# Name <- 'H3122-LR_ctrl-1'

setwd("/disk0/sm/bulk/test/output")
cfiles <- list.files("/disk0/sm/bulk/test/output", pattern='ReadsPerGene.out.tab', recursive=TRUE)
cfiles

# non strand option. if strand specific, columns = c(1,4)
counts <- readDGE(cfiles,columns=c(1,4))
a <- unlist(strsplit(colnames(counts), split='/'))
colnames(counts) <- a[seq(1,length(a),by=2)]

noint = rownames(counts) %in% 
c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")
# CRITICAL STEP In edgeR, it is recommended to remove features without at least 1 read per million in n of the 
# samples, where n is the size of the smallest group of replicates (here, n = 3 for the knockdown group).
cpms = cpm(counts)
keep = rowSums(cpms > 1) >= 1 & !noint
counts = counts[keep,]

######################################### factor 수동설정

Tissue <- factor(c("control","control","control","case","case","case"))
design <- model.matrix(~0+Tissue)
y <- DGEList(counts = counts, group = Tissue)
y <- calcNormFactors(y)
normalized_cpm = cpm(y)
normalized_cpm_log2 = log2(cpm(y)+1)

write.table(file="/disk0/sm/bulk/test/output/cpm_normalized_result.txt",normalized_cpm,sep='\t',col.names=NA,quote=FALSE)
write.table(file="/disk0/sm/bulk/test/output/cpm_normalized_log2_result.txt",normalized_cpm_log2,sep='\t',col.names=NA,quote=FALSE)

pdf("/disk0/sm/bulk/test/output/MDSplot.pdf")
plotMDS(y)
dev.off()

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=c(-1,+1))
topTags(lrt)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(file="/disk0/sm/bulk/test/output/edgeR_CDA_vs_LR_Total.txt",out,sep='\t',col.names=NA,quote=FALSE)

keep2 <- out$table$FDR <= 0.05 & abs(out$table$logFC) >= 1
write.table(file="/disk0/sm/bulk/test/output/edgeR_CDA_vs_LR_DEG.txt",out[keep2,],sep='\t',col.names=NA,quote=FALSE)

write.table(file="/disk0/sm/bulk/test/output/edgeR_DEG_expression.txt",normalized_cpm_log2[rownames(out[keep2,]),],sep='\t',col.names=NA,quote=FALSE)
