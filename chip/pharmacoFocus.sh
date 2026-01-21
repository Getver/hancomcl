# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KORV2 KORV2_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KORV1 KORV1_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 PMDA PMDA_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 PangenomiX PangenomiX_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KBA_B 0926

#!/bin/bash

# sudo chmod -R 775 /disk2/hancomcl/PMDA/

## library PMDA 적용(전에 하던대로 했더니 돌아감)

KEY=$1
LIBRARY=$2

mkdir /workdir/${KEY}
cat /code_SM/header.txt /workdir/${KEY}.txt | sort | uniq > /workdir/${KEY}/${KEY}_celfile.txt
rm /workdir/${KEY}.txt


library_folder='/library_QC/library/Axiom_PharmacoFocus.r6'

CDF_file=${library_folder}'/Axiom_PharmacoFocus.r6.cdf'
QCC_file=${library_folder}'/Axiom_PharmacoFocus.r6.qcc'
QCA_file=${library_folder}'/Axiom_PharmacoFocus.r6.qca'

ARG_file=${library_folder}'/Axiom_PharmacoFocus_96orMore_Step2.r6.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml'
CNARG_file=${library_folder}'/Axiom_PharmacoFocus.r6.apt-genotype-axiom.AxiomCN_PS1.apt2.xml'
CNVmix_file=${library_folder}'/Axiom_PharmacoFocus.r6.apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml'
CN_models_file=${library_folder}'/Axiom_PharmacoFocus.r6.cn_models'

SPECIALSNP_file=${library_folder}'/Axiom_PharmacoFocus.r6.specialSNPs'
PS2SNP_file=${library_folder}'/Axiom_PharmacoFocus.r6.ps2multisnp_map.ps'

TRANSLATION_file=${library_folder}'/Axiom_PharmacoFocus.r6.translation'
METABOLIZER_file=${library_folder}'/Axiom_PharmacoFocus.r6.metabolizer'
DCANNOT_file=${library_folder}'/Axiom_PharmacoFocus.r6.dc_annot.csv'
DB_file=${library_folder}'/Axiom_PharmacoFocus.na36.r6.a2.annot.db'

## 1. DQC

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-geno-qc \
--cdf-file ${CDF_file} \
--qcc-file ${QCC_file} \
--qca-file ${QCA_file} \
--out-dir /workdir/${KEY}/01_Genotype \
--out-file /workdir/${KEY}/01_Genotype/apt-geno-qc.txt \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt

echo "["${LIBRARY}"]"${KEY}" - DQC finish! FOLDER : 01_Genotype"


## 2. Generate summary signals for all probesets

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/Axiom_PharmacoFocus.r6 \
--arg-file ${CNARG_file} \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--out-dir /workdir/${KEY}/03_CN_input \
--log-file /workdir/${KEY}/03_CN_input/apt2-axiom.log \
--allele-summaries-file /workdir/${KEY}/03_CN_input/AxiomGT1.summary.a5 
# --summary-a5-output False \

echo "["${LIBRARY}"]"${KEY}" - CN input finish! FOLDER : 03_CN_input"


## 3. Run the copy number analysis

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-copynumber-axiom-cnvmix \
--analysis-files-path ${library_folder} \
--arg-file ${CNVmix_file} \
--reference-file ${CN_models_file} \
--mapd-max 0.35 \
--waviness-sd-max 0.1 \
--summary-file /workdir/${KEY}/03_CN_input/AxiomGT1.summary.a5 \
--report-file /workdir/${KEY}/03_CN_input/AxiomGT1.report.txt \
--out-dir /workdir/${KEY}/04_CN \
--log-file /workdir/${KEY}/04_CN/apt-copynumber-axiom.log

echo "["${LIBRARY}"]"${KEY}" - CN finish! FOLDER : 04_CN"


## 4. final genotyping

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path ${library_folder} \
--arg-file ${ARG_file} \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--out-dir /workdir/${KEY}/01_Genotype \
--batch-folder /workdir/${KEY}/01_Genotype \
--genotyping-node:snp-posteriors-output true \
--multi-genotyping-node:multi-posteriors-output true \
--allele-summaries true \
--do-rare-het-adjustment true \
--multi-posteriors-output true \
--log-file /workdir/${KEY}/01_Genotype/apt2-axiom.log

echo "["${LIBRARY}"]"${KEY}" - Final genotyping finish! FOLDER : 01_Genotype"


## 5. Generate SNP QC metrics

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-metrics \
--posterior-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.txt \
--batch-folder /workdir/${KEY}/01_Genotype \
--report-file /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt \
--special-snps ${SPECIALSNP_file} \
--y-restrict 0.2 \
--min-genotype-freq-samples 20 \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--log-file /workdir/${KEY}/02_SNPolisher/ps_metrics.log 
# --use-multi-allele true \
# --multi-posterior-file /disk2/kb/sm/PMDA/${name}/genotypes/AxiomGT1.snp-posteriors.multi.txt \
# --multi-metrics-file /disk2/kb/sm/PMDA/${name}/SNPolisher/metrics.multi.txt \

echo "["${LIBRARY}"]"${KEY}" - SNP QC metrics finish! FOLDER : 02_SNPolishers"


## 6. Run SNP classification

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-classification \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--ps2snp-file ${PS2SNP_file} \
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
--output-dir /workdir/${KEY}/02_SNPolisher \
--log-file /workdir/${KEY}/02_SNPolisher/ps_classification.log \

echo "["${LIBRARY}"]"${KEY}" - SNP classification finish! FOLDER : 02_SNPolishers"


## 7. STAR allele translation


/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt2-pgx-translation \
--batch-folder /workdir/${KEY}/01_Genotype \
--cn-region-calls-file /workdir/${KEY}/04_CN/AxiomCNVMix.cnregioncalls.txt \
--marker-list-file /workdir/${KEY}/02_SNPolisher/Recommended.ps \
--analysis-files-path ${library_folder} \
--translate-file ${TRANSLATION_file} \
--metabolizer-file ${METABOLIZER_file} \
--annotation-file ${DCANNOT_file} \
--out-dir /workdir/${KEY}/05_Translation \
--base-name-prefix ${KEY} \
--log-file /workdir/${KEY}/05_Translation/${KEY}.log \
--use-first-dup-allele-def true 

echo "["${LIBRARY}"]"${KEY}" - STAR allele translation finish! FOLDER : 05_Translation"


/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-format-result \
--calls-file /workdir/${KEY}/01_Genotype/AxiomGT1.calls.txt \
--annotation-file ${DB_file} \
--export-vcf-file /workdir/${KEY}/06_vcf/${KEY}.vcf \
--log-file /workdir/${KEY}/06_vcf/apt-format-result.log \
--snp-identifier-column Affy_SNP_ID ## 바꿀 수 있음


echo "["${LIBRARY}"]"${KEY}" - final VCF finish! FOLDER : 06_vcf"

bgzip /workdir/${KEY}/06_vcf/${KEY}.vcf
tabix -p vcf /workdir/${KEY}/06_vcf/${KEY}.vcf.gz

grep -v '^#' /workdir/${KEY}/01_Genotype/apt-geno-qc.txt | awk '{print $1"\t"$18}' > /workdir/${KEY}/01_Genotype/DQC_rate.txt
grep -v '^#' /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$19}' > /workdir/${KEY}/01_Genotype/QC_rate.txt
paste /workdir/${KEY}/01_Genotype/QC_rate.txt /workdir/${KEY}/01_Genotype/DQC_rate.txt | awk '{$6=""; print $0}'> /workdir/${KEY}/01_Genotype/QC_table.txt

python3 /code_SM/make_table.py ${KEY} ${LIBRARY}

chmod -R 777 /workdir/${KEY}
