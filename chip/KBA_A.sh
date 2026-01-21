# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KORV2 KORV2_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KORV1 KORV1_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 PMDA PMDA_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 PangenomiX PangenomiX_test
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl:/disk2/hancomcl -v /disk0/sm/chip/docker/library:/library_QC/library sm:chip_0.2 KBA_B 0926

#!/bin/bash
# sudo chmod -R 775 /disk2/hancomcl/KORV2/

## library Axiom_KBA_v2_A 적용

KEY=$1
LIBRARY=$2
FOLDER='Axiom_KBA_v2_A'

mkdir /workdir/${KEY}
cat /code_SM/header.txt /workdir/${KEY}.txt | sort | uniq > /workdir/${KEY}/${KEY}_celfile.txt
rm /workdir/${KEY}.txt


## DQC

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-geno-qc \
--cdf-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.cdf \
--qcc-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.qcc \
--qca-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.qca \
--out-dir /workdir/${KEY}/01_Genotype \
--out-file /workdir/${KEY}/01_Genotype/apt-geno-qc.txt \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt

echo "["${LIBRARY}"]"${KEY}" - DQC finish! FOLDER : 01_Genotype"


## genotype calling

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/${FOLDER} \
--arg-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
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

echo "["${LIBRARY}"]"${KEY}" - SNP QC metrics finish! FOLDER : 02_SNPolisher"


## Classification

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/ps-classification \
--species-type human \
--metrics-file /workdir/${KEY}/02_SNPolisher/metrics.txt \
--output-dir /workdir/${KEY}/02_SNPolisher \
--log-file /workdir/${KEY}/02_SNPolisher/ps-classification.log \
--multi-metrics-file /workdir/${KEY}/02_SNPolisher/multi-metrics.txt

echo "["${LIBRARY}"]"${KEY}" - SNP classification finish! FOLDER : 02_SNPolisher"


## SSP

/library_QC/tool/simple-ssps \
--posterior-file /workdir/${KEY}/01_Genotype/AxiomGT1.snp-posteriors.txt \
--performance-file /workdir/${KEY}/02_SNPolisher/Ps.performance.txt  \
--output-dir /workdir/${KEY}/01_Genotype \
--log-file /workdir/${KEY}/01_Genotype/ssplog.txt 

echo "["${LIBRARY}"]"${KEY}" - SSP finish! FOLDER : 01_Genotype"


## Re-genotyping using SSP

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /library_QC/library/${FOLDER}/ \
--arg-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
--out-dir /workdir/${KEY}/04_Regenotype \
--cel-files /workdir/${KEY}/${KEY}_celfile.txt \
--genotyping-node:snp-priors-input-file /workdir/${KEY}/01_Genotype/AxiomGT1.*.models \
--log-file /workdir/${KEY}/04_Regenotype/apt-regenotype-axiom.log \
--do-rare-het-adjustment True \
--artifact-reduction-output-trustcheck True \
--genotyping-node:snp-posteriors-output True \
--allele-summaries \
--allele-summaries-file /workdir/${KEY}/04_Regenotype/AxiomGT1.summary.txt \
--summary-a5-output False \
--multi-posteriors-output True \
--multi-posteriors-output-file /workdir/${KEY}/04_Regenotype/AxiomGT1.snp-posteriors.multi.txt

echo "["${LIBRARY}"]"${KEY}" - Re-genotyping finish! FOLDER : 04_Regenotype"


## Advanced normalization

cd /library_QC/advnorm-1.1

bash advnorm.sh \
--summary-file /workdir/${KEY}/04_Regenotype/AxiomGT1.summary.txt \
--calls-file /workdir/${KEY}/04_Regenotype/AxiomGT1.calls.txt \
--report-file /workdir/${KEY}/04_Regenotype/AxiomGT1.report.txt \
--trustcheck-file /workdir/${KEY}/04_Regenotype/AxiomGT1.trustcheck.txt \
--analysis-files-path /library_QC/library/${FOLDER} \
--snp-priors-file /workdir/${KEY}/01_Genotype/AxiomGT1.*.models \
--special-snps-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.specialSNPs \
--ps2snp-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.ps2snp_map.ps \
--output-dir /workdir/${KEY}/05_Normalization \
--snp-specific-param-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.r1.probeset_genotyping_parameters.txt

echo "["${LIBRARY}"]"${KEY}" - Advanced normalization finish! FOLDER : 05_Normalization"


## vcf로 변환

/library_QC/tool/apt_2.12.0_linux_64_x86_binaries/bin/apt-format-result \
--calls-file /workdir/${KEY}/01_Genotype/*.calls.txt \
--annotation-file /library_QC/library/${FOLDER}/Axiom_KBA_v2_A.na36.r1.a1.annot.db \
--export-vcf-file /workdir/${KEY}/03_vcf/${KEY}.vcf \
--log-file /workdir/${KEY}/03_vcf/apt-format-result.log

echo "["${LIBRARY}"]"${KEY}" - VCF finish! FOLDER : 03_vcf"


bgzip /workdir/${KEY}/03_vcf/${KEY}.vcf
tabix -p vcf /workdir/${KEY}/03_vcf/${KEY}.vcf.gz


grep -v '^[#]' /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt | awk '{ print $1"\t"$3 }' > /workdir/${KEY}/CallRate.txt
grep -v '^[#]' /workdir/${KEY}/04_Regenotype/AxiomGT1.report.txt | awk '{ print $1"\t"$3 }' > /workdir/${KEY}/CallRate_re.txt
grep -v '^[#]' /workdir/${KEY}/05_Normalization/AxiomGT1.report.txt | awk '{ print $1"\t"$3 }' > /workdir/${KEY}/CallRate_normal.txt
grep -v '^[#]' /workdir/${KEY}/01_Genotype/AxiomGT1.report.txt | awk '{ print $1"\t"$2 }' > /workdir/${KEY}/Gender.txt
grep -v '^[#]' /workdir/${KEY}/01_Genotype/apt-geno-qc.txt | awk '{ print $1"\t"$18 }' > /workdir/${KEY}/DQC_Rate.txt

python3 /code_SM/make_table.py ${KEY} ${LIBRARY}

chmod -R 777 /workdir/${KEY}

# platenum=$(less cel_list2.txt | awk -F '_' '{ print $3 }' | uniq | wc -l)

# Rscript ../plot.R $1 $platenum

# cd /disk2/kb/sm/chip/QC

# python3 QC_table.py $1

# python3 Recommend.py $1


# sudo chmod 775 /disk2/hancomcl/KORV2/*1204*.CEL
# ls /disk2/hancomcl/JPMI/*.CEL | awk -F '_' '{ print $3 }' | uniq -c
# ls /disk2/hancomcl/KORV1/2023/*.CEL | awk -F '_' '{ print $2 }' | uniq -c
# ls /disk2/hancomcl/KORV2/*.CEL | awk -F '_' '{ print $3 }' | uniq -c
# ls /disk2/hancomcl/PMDA/*.CEL | awk -F '_' '{ print $3 }' | uniq -c
# ll --time-style full-iso /disk2/hancomcl/KORV2/*.CEL > test.sumin
# less test.sumin | awk -F '[ _]' '{ print $6"\t"$11 }' | sort -k2 | uniq
# ls /disk2/hancomcl/KORV2/*.CEL | shuf -n 2000 | sort -t_ -k3 > test5/test5.txt


# cat test11.txt test22.txt | sort | uniq -u > test1.txt
