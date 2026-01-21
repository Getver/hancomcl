# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KORV2 KORV2_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KORV1 KORV1_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 PMDA PMDA_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 PangenomiX PangenomiX_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/library:/library_QC/library sm:chip_0.2 KBA_B 0926


#!/bin/bash
# 647359


# 50K chip
## library r1 적용

KEY=$1
LIBRARY=$2

mkdir /workdir/${KEY}
cat /code_SM/header.txt /workdir/${KEY}.txt | sort | uniq > /workdir/${KEY}/${KEY}_celfile.txt
rm /workdir/${KEY}.txt


## DQC

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-geno-qc \
--cdf-file /library_QC/library/Axiom_HCLV1/Axiom_HCLV1.r1.cdf \
--qcc-file /library_QC/library/Axiom_HCLV1/Axiom_HCLV1.r1.qcc \
--qca-file /library_QC/library/Axiom_HCLV1/Axiom_HCLV1.r1.qca \
--out-dir /workdir/${KEY}/01_Genotype \
--out-file /workdir/${KEY}/01_Genotype/apt-geno-qc.txt \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt

echo "["${LIBRARY}"]"${KEY}" - DQC finish! FOLDER : 01_Genotype"


## genotype calling

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/Axiom_HCLV1 \
--arg-file /library_QC/library/Axiom_HCLV1/Axiom_HCLV1_96orMore_Step2.r1.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
--out-dir /workdir/${KEY}/01_Genotype \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--log-file /workdir/${KEY}/01_Genotype/apt-genotype-axiom.log \
--do-rare-het-adjustment True \
--artifact-reduction-output-trustcheck False \
--genotyping-node:snp-posteriors-output True \
--allele-summaries \
--allele-summaries-file /workdir/${KEY}/01_Genotype/AxiomGT1.summary.txt \
--summary-a5-output False 

echo "["${LIBRARY}"]"${KEY}" - Genotype calling finish!! FOLDER : 01_Genotype"


## QC matrix

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-metrics \
--posterior-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.txt \
--call-file /workdir/${KEY}/01_Genotype/AxiomGT1.calls.txt \
--summary-file /workdir/${KEY}/01_Genotype/AxiomGT1.summary.txt \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--log-file /workdir/${KEY}/02_SNPolisher/ps-metrics.log 

echo "["${LIBRARY}"]"${KEY}" - QC matrix finish! FOLDER : 02_SNPolisher"


## Classification

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-classification \
--species-type human \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--output-dir /workdir/${KEY}/02_SNPolisher \
--log-file /workdir/${KEY}/02_SNPolisher/ps-classification.log 

echo "["${LIBRARY}"]"${KEY}" - Classification finish! FOLDER : 02_SNPolisher"


## vcf로 변환

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-format-result \
--calls-file /workdir/${KEY}/01_Genotype/*.calls.txt \
--annotation-file /library_QC/library/Axiom_HCLV1/Annotations_hg38/Axiom_HCLV1.na36.r1.a1.annot.db \
--export-vcf-file /workdir/${KEY}/03_vcf/${KEY}.vcf \
--log-file /workdir/${KEY}/03_vcf/apt-format-result.log

echo "["${LIBRARY}"]"${KEY}" - VCF finish! FOLDER : 03_vcf"


bgzip /workdir/${KEY}/03_vcf/${KEY}.vcf
tabix -p vcf /workdir/${KEY}/03_vcf/${KEY}.vcf.gz


grep -v '^#' /workdir/${KEY}/01_Genotype/apt-geno-qc.txt | awk '{print $1"\t"$18}' > /workdir/${KEY}/01_Genotype/DQC_rate.txt
grep -v '^#' /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$19}' > /workdir/${KEY}/01_Genotype/QC_rate.txt
paste /workdir/${KEY}/01_Genotype/QC_rate.txt /workdir/${KEY}/01_Genotype/DQC_rate.txt | awk '{$6=""; print $0}'> /workdir/${KEY}/01_Genotype/QC_table.txt

python3 /code_SM/make_table.py ${KEY} ${LIBRARY}

chmod -R 777 /workdir/${KEY}

