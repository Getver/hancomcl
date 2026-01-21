library(ggplot2)
library(reshape2)
library(SNPolisher)
library(R6)

args<-commandArgs(trailingOnly=TRUE);
KEY <- args[1];
LIBRARY <- args[2];

# if (LIBRARY == 'KORV1' || LIBRARY == 'KORV2'){
#     Ps_Visualization (
#         pidFile=paste0('/workdir/',KEY,'/02_SNPolisher/Recommended.ps'), 
#         summaryFile=paste0('/workdir/',KEY,'/01_Genotype/AxiomGT1.summary.txt'), 
#         callFile=paste0('/workdir/',KEY,'/01_Genotype/AxiomGT1.calls.txt'), 
#         confidenceFile=paste0('/workdir/',KEY,'/01_Genotype/AxiomGT1.confidences.txt'), 
#         posteriorFile=paste0('/workdir/',KEY,'/01_Genotype/AxiomGT1.snp-posteriors.txt'), 
#         output.File=paste0('/workdir/',KEY,'/ClusteringPlot_6.pdf'),
#         refFile=NULL,
#         plot.prior=TRUE,
#         priorFile=NULL,
#         sampleFile=NULL,
#         max.num.SNP.draw=6
#     )
# } else 

if (LIBRARY == 'PMDA'){
    CN_Visualization(
        CNregioncallsFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cnregioncalls.txt'),
        CNpriorFile = '/library_QC/library/Axiom_PMDA_Plus.Analysis_r7.1/Axiom_PMDA.r7.with_mPCR.cn_priors',
        CNposteriorFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cn_posteriors.txt'), 
        CNdetailsFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cnregions.details.txt'),
        CNtruthFile = NULL, 
        output.dir=paste0('/workdir/', KEY),
        output.File = 'CNV_plot.pdf', 
        plot.plate.effects = FALSE
    )
} else if (LIBRARY == 'PangenomiX'){
    CN_Visualization(
        CNregioncallsFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cnregioncalls.txt'),
        CNpriorFile = '/library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.cn_priors',
        CNposteriorFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cn_posteriors.txt'), 
        CNdetailsFile = paste0('/workdir/', KEY, '/05_CN/AxiomCNVMix.cnregions.details.txt'),
        CNtruthFile = NULL, 
        output.dir=paste0('/workdir/', KEY),
        output.File = 'CNV_plot.pdf', 
        plot.plate.effects = FALSE
    )
} else if (LIBRARY == 'pharmacoFocus'){
    CN_Visualization(
        CNregioncallsFile = paste0('/workdir/', KEY, '/04_CN/AxiomCNVMix.cnregioncalls.txt'),
        CNpriorFile = '/library_QC/library/Axiom_PharmacoFocus.r6/Axiom_PharmacoFocus.r6.cn_priors',
        CNposteriorFile = paste0('/workdir/', KEY, '/04_CN/AxiomCNVMix.cn_posteriors.txt'), 
        CNdetailsFile = paste0('/workdir/', KEY, '/04_CN/AxiomCNVMix.cnregions.details.txt'),
        CNtruthFile = NULL, 
        output.dir=paste0('/workdir/', KEY),
        output.File = 'CNV_plot.pdf', 
        plot.plate.effects = FALSE
    )
} 
