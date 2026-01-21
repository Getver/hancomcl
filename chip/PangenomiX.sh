# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KORV2 KORV2_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KORV1 KORV1_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 PMDA PMDA_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 PangenomiX PangenomiX_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KBA_B 0926


#!/bin/bash

# sudo chmod -R 775 /disk2/hancomcl/PMDA/

## library PMDA 적용(전에 하던대로 했더니 돌아감)

KEY=$1
LIBRARY=$2

mkdir /workdir/${KEY}
cat /code_SM/header.txt /workdir/${KEY}.txt | sort | uniq > /workdir/${KEY}/${KEY}_celfile.txt
rm /workdir/${KEY}.txt


## DQC

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-geno-qc \
--cdf-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.cdf \
--qcc-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.qcc \
--qca-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.qca \
--out-dir /workdir/${KEY}/01_Genotype \
--out-file /workdir/${KEY}/01_Genotype/apt-geno-qc.txt \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt

echo "["${LIBRARY}"]"${KEY}" - DQC finish! FOLDER : 01_Genotype"


## genotype calling
# <Parameter name="artifact-reduction-clip-pvcam" analysis="artifact-reduction-node" currentValue="0.43" />
/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus \
--arg-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX_96orMore_Step2.r1.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
--out-dir /workdir/${KEY}/01_Genotype \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--log-file /workdir/${KEY}/01_Genotype/apt-genotype-axiom.log \
--do-rare-het-adjustment True \
--artifact-reduction-output-trustcheck True \
--genotyping-node:snp-posteriors-output True \
--allele-summaries \
--allele-summaries-file /workdir/${KEY}/01_Genotype/AxiomGT1.summary.txt \
--summary-a5-output False \
--multi-posteriors-output True \
--multi-posteriors-output-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.multi.txt

echo "["${LIBRARY}"]"${KEY}" - Genotype calling finish! FOLDER : 01_Genotype"


## QC matrix

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-metrics \
--posterior-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.txt \
--call-file /workdir/${KEY}/01_Genotype/AxiomGT1.calls.txt \
--summary-file /workdir/${KEY}/01_Genotype/AxiomGT1.summary.txt \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--log-file /workdir/${KEY}/02_SNPolisher/ps-metrics.log \
--multi-posterior-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.multi.txt

echo "["${LIBRARY}"]"${KEY}" - QC matrix finish! FOLDER : 02_SNPolisher"


## Classification

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-classification \
--species-type human \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--output-dir /workdir/${KEY}/02_SNPolisher/ \
--log-file /workdir/${KEY}/02_SNPolisher/ps-classification.log \
--multi-metrics-file /workdir/${KEY}/02_SNPolisher/multi-metrics.txt

echo "["${LIBRARY}"]"${KEY}" - Classification finish! FOLDER : 02_SNPolisher"


## vcf로 변환

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-format-result \
--calls-file /workdir/${KEY}/01_Genotype/AxiomGT1.calls.txt \
--annotation-file /library_QC/library/Axiom_PangenomiX.r1/Axiom_PangenomiX.na36.r1.a1.annot.db \
--export-vcf-file /workdir/${KEY}/03_vcf/${KEY}.vcf \
--log-file /workdir/${KEY}/03_vcf/apt-format-result.log

echo "["${LIBRARY}"]"${KEY}" - VCF finish! FOLDER : 03_vcf"


## 3. Generate summary signals for all probesets

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus \
--arg-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.apt-genotype-axiom.AxiomCN_PS1.apt2.xml \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--out-dir /workdir/${KEY}/04_CN_input \
--log-file /workdir/${KEY}/04_CN_input/apt2-axiom.log \
--allele-summaries-file /workdir/${KEY}/04_CN_input/AxiomGT1.summary.a5 
# --summary-a5-output False \

echo "["${LIBRARY}"]"${KEY}" - CN input finish! FOLDER : 04_CN_input"


## 4. Run the copy number analysis

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-copynumber-axiom-cnvmix \
--analysis-files-path /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus \
--arg-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml \
--reference-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.cn_models \
--mapd-max 0.35 \
--waviness-sd-max 0.1 \
--summary-file /workdir/${KEY}/04_CN_input/AxiomGT1.summary.a5 \
--report-file /workdir/${KEY}/04_CN_input/AxiomGT1.report.txt \
--out-dir /workdir/${KEY}/05_CN \
--log-file /workdir/${KEY}/05_CN/apt-copynumber-axiom.log

echo "["${LIBRARY}"]"${KEY}" - CN finish! FOLDER : 05_CN"


## 5. final genotyping

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--copynumber-probeset-calls-file /workdir/${KEY}/05_CN/AxiomCNVMix.cnpscalls.txt \
--analysis-files-path /library_QC/library/Axiom_PangenomiX.r1 \
--arg-file /library_QC/library/Axiom_PangenomiX.r1/Axiom_PangenomiX_96orMore_Step2.r1.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--out-dir /workdir/${KEY}/06_Genotype \
--batch-folder /workdir/${KEY}/06_Genotype \
--genotyping-node:snp-posteriors-output true \
--multi-genotyping-node:multi-posteriors-output true \
--allele-summaries true \
--do-rare-het-adjustment true \
--multi-posteriors-output true \
--log-file /workdir/${KEY}/06_Genotype/apt2-axiom.log

echo "["${LIBRARY}"]"${KEY}" - Final genotyping finish! FOLDER : 06_Genotype"


## 6. Generate SNP QC metrics

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-metrics \
--posterior-file /workdir/${KEY}/06_Genotype/AxiomGT1.snp-posteriors.txt \
--batch-folder /workdir/${KEY}/06_Genotype \
--report-file /workdir/${KEY}/06_Genotype/AxiomGT1.report.txt \
--special-snps /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.specialSNPs \
--y-restrict 0.2 \
--min-genotype-freq-samples 20 \
--metrics-file /workdir/${KEY}/07_SNPolisher/metrics.txt \
--log-file /workdir/${KEY}/07_SNPolisher/ps_metrics.log 
# --use-multi-allele true \
# --multi-posterior-file /disk2/kb/sm/PMDA/${name}/genotypes/AxiomGT1.snp-posteriors.multi.txt \
# --multi-metrics-file /disk2/kb/sm/PMDA/${name}/SNPolisher/metrics.multi.txt \

echo "["${LIBRARY}"]"${KEY}" - SNP QC metrics finish! FOLDER : 07_SNPolishers"


## 7. Run SNP classification

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-classification \
--metrics-file /workdir/${KEY}/07_SNPolisher/metrics.txt \
--ps2snp-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.ps2multisnp_map.ps \
--species-type Human \
--cr-cutoff 95 \
--fld-cutoff 3.6 \
--het-so-cutoff -0.1 \
--het-so-XChr-cutoff -0.1 \
--het-so-otv-cutoff -0.3 \
--hom-ro-1-cutoff 0.6 \
--hom-ro-2-cutoff 0.3 \
--hom-ro-3-cutoff -0.9 \
--hom-ro true \
--num-minor-allele-cutoff 2 \
--hom-ro-hap-1-cutoff 0.1 \
--hom-ro-hap-1-XChr-cutoff 0.1 \
--hom-ro-hap-1-MTChr-cutoff 0.4 \
--hom-ro-hap-2-cutoff -0.9 \
--hom-ro-hap-2-XChr-cutoff 0.05 \
--hom-ro-hap-2-MTChr-cutoff 0.2 \
--hom-hap-X-cutoff -1 \
--hom-hap-Y-lower-cutoff 1 \
--hom-hap-Y-upper-cutoff 1 \
--CN0-hap-X-cutoff -1 \
--CN0-hap-Y-cutoff -1 \
--CN0-dip-X-cutoff -1 \
--CN0-dip-Y-cutoff -1 \
--aaf-XChr-cut 0.36 \
--fld-XChr-cut 4 \
--homfld-XChr-cut 6.5 \
--homfld-YChr-cut 6.5 \
--min-YChr-samples-cut 5 \
--priority-order PolyHighResolution,NoMinorHom,MonoHighResolution,OTV,UnexpectedGenotypeFreq,CallRateBelowThreshold,Other,OtherMA \
--recommended PolyHighResolution,NoMinorHom,MonoHighResolution,Hemizygous \
--genotype-p-value-cutoff 0.000001 \
--hom-mma-cutoff 8.5 \
--fld-ma-cutoff 4.2 \
--fld-ma-2-cutoff 4.2 \
--min-fld-ma-cutoff 3.2 \
--min-fld-ma-2-cutoff 3.2 \
--het-so-ma-2-cutoff -0.2 \
--hom-ro-ma-cutoff 0.2 \
--hom-ro-ma-2-cutoff 0.3 \
--hom-ro-ma-1-cutoff 0.3 \
--priority-order-multi-allele PolyHighResolution,NoMinorHom,MonoHighResolution,Hemizygous,UnexpectedGenotypeFreq,CallRateBelowThreshold,Other,OtherMA \
--best-cr-ma-cutoff 90.0 \
--output-dir /workdir/${KEY}/07_SNPolisher \
--log-file /workdir/${KEY}/07_SNPolisher/ps_classification.log \

echo "["${LIBRARY}"]"${KEY}" - SNP classification finish! FOLDER : 07_SNPolishers"


## STAR allele translation

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt2-pgx-translation \
--batch-folder /workdir/${KEY}/06_Genotype \
--cn-region-calls-file /workdir/${KEY}/05_CN/AxiomCNVMix.cnregioncalls.txt \
--marker-list-file /workdir/${KEY}/07_SNPolisher/Recommended.ps \
--analysis-files-path /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus \
--translate-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.translation \
--metabolizer-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.metabolizer \
--annotation-file /library_QC/library/Axiom_PangenomiX_Analysis.r1_plus/Axiom_PangenomiX.r1.dc_annot.csv \
--out-dir /workdir/${KEY}/08_Translation \
--base-name-prefix ${KEY} \
--log-file /workdir/${KEY}/08_Translation/${KEY}.log \
--use-first-dup-allele-def true 

echo "["${LIBRARY}"]"${KEY}" - STAR allele translation finish! FOLDER : 08_Translation"


## generate final VCF files

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-format-result \
--batch-folder /workdir/${KEY}/06_Genotype \
--cn-region-calls-file /workdir/${KEY}/05_CN/AxiomCNVMix.cnregioncalls.txt \
--snp-list-file /workdir/${KEY}/02_SNPolisher/Recommended.ps \
--annotation-file /library_QC/library/Axiom_PangenomiX.r1/Axiom_PangenomiX.na36.r1.a1.annot.db \
--snp-identifier-column custom_rsid \
--export-chr-shortname true \
--export-vcf-file /workdir/${KEY}/09_Export/${KEY}.vcf \
--log-file /workdir/${KEY}/09_Export/${KEY}.log

echo "["${LIBRARY}"]"${KEY}" - final VCF finish! FOLDER : 09_Export"


bgzip /workdir/${KEY}/03_vcf/${KEY}.vcf
tabix -p vcf /workdir/${KEY}/03_vcf/${KEY}.vcf.gz
bgzip /workdir/${KEY}/09_Export/${KEY}.vcf
tabix -p vcf /workdir/${KEY}/09_Export/${KEY}.vcf.gz

grep -v '^#' /workdir/${KEY}/01_Genotype/apt-geno-qc.txt | awk '{print $1"\t"$18}' > /workdir/${KEY}/01_Genotype/DQC_rate.txt
grep -v '^#' /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$19}' > /workdir/${KEY}/01_Genotype/QC_rate.txt
paste /workdir/${KEY}/01_Genotype/QC_rate.txt /workdir/${KEY}/01_Genotype/DQC_rate.txt | awk '{$6=""; print $0}'> /workdir/${KEY}/01_Genotype/QC_table.txt

python3 /code_SM/make_table.py ${KEY} ${LIBRARY}

chmod -R 777 /workdir/${KEY}

