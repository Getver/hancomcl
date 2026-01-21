# docker run -it \
# -v /disk2/kb/sm/bulk/fastq/test:/data \
# -v /disk2/kb/exome_sequencing:/knownSite \
# -v /disk2/references/hg38/UCSC/DNA/all_chromosomes:/reference \
# -v /disk2/kb/sm/T7/test:/root \
# sm:0.0.1

OUT_DIR=
INPUT_DIR=
F1=
F2=
REFERENCE_FILE=
KEY=
RG="'@RG\\tID:${KEY}\\tSM:${KEY}\\tPL:DNBSEQ\\tPU:flowcell_barcode\'"


fastqc -o ${OUT_DIR} ${F1} ${F2}


## bam 만들기
bwa mem -t 32 ${REFERENCE_FILE} ${F1} ${F2} -R ${RG} | samtools view -h -b -o ${KEY}.bam

## sortsam 써서 sort_bam 만들기
gatk --java-options "-Xmx4G" SortSam \
    -I ${KEY}.bam \
    -O ${KEY}_sorted.bam -SO coordinate \
    --VALIDATION_STRINGENCY SILENT \
    --MAX_RECORDS_IN_RAM 5000000

samtools index ${KEY}_sorted.bam

## on-target 얻고싶을때 쓰기(즉 WES일때만 사용하기)
# java -jar /disk2/kb/sm/cancer/tool/picard.jar CollectHsMetrics \
#     --INPUT ${KEY}_sorted.bam \
#     --OUTPUT ${KEY}_hsmetrics.txt \
#     --BAIT_INTERVALS /disk2/kb/G400/target.interval_list \
#     --TARGET_INTERVALS /disk2/kb/G400/target.interval_list

## qualimap(sorted 하고나서 해야함)
# qualimap --java-mem-size=80G bamqc \
#     -bam ${KEY}_sorted.bam \
#     --outfile ${KEY}

## markduplicates
gatk --java-options "-Xmx4G" MarkDuplicates \
    -I ${KEY}_sorted.bam \
    -O ${KEY}_dup.bam \
    -M ${KEY}_dup.metrics \
    --VALIDATION_STRINGENCY SILENT \
    --MAX_RECORDS_IN_RAM 5000000

samtools index ${KEY}_dup.bam

## BaseRecalibrator
gatk --java-options "-Xmx4G" BaseRecalibrator -R ${REFERENCE_FILE} -I ${KEY}_dup.bam \
    --known-sites /knownSite/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /knownSite/hapmap_3.3.hg38.vcf.gz \
    --known-sites /knownSite/1000G_omni2.5.hg38.vcf.gz \
    --known-sites /knownSite/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites /knownSite/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites /knownSite/Homo_sapiens.dbsnp155_final.vcf.gz \
    -O ${KEY}_dup.recal

gatk --java-options "-Xmx4G" ApplyBQSR -R ${REFERENCE_FILE} -I ${KEY}_dup.bam --bqsr-recal-file ${KEY}_dup.recal -O ${KEY}_recal.bam \
    --create-output-bam-index false

samtools index ${KEY}_recal.bam

## multiplemetrics 얻고싶을때 쓰기
# gatk --java-options "-Xmx4G" CollectMultipleMetrics \
#     -R ${REFERENCE_FILE} \
#     -I ${KEY}_recal.bam \
#     -O ${KEY}_aligned \
#     -AS true \
#     --PROGRAM "CollectAlignmentSummaryMetrics" \
#     --PROGRAM "CollectInsertSizeMetrics" \
#     --PROGRAM "CollectSequencingArtifactMetrics" \
#     --PROGRAM "CollectGcBiasMetrics" \
#     --PROGRAM "QualityScoreDistribution"

## haplotypecaller(population vcf를 만들고싶으면 gvcf로 만들고 그냥 한명 vcf로 만들고싶으면 -ERC 옵션을 빼면 됨)
gatk --java-options "-Xms8G" HaplotypeCaller \
    -R ${REFERENCE_FILE} \
    -I ${KEY}_recal.bam \
    -O ${OUT_DIR}/gvcf/${KEY}.gvcf.gz \
    -D /knownSite/Homo_sapiens.dbsnp155_final.vcf.gz \
    -ERC GVCF

tabix -p vcf ${OUT_DIR}/gvcf/${KEY}.gvcf.gz




## 여기부턴 population vcf를 만들기 위한 과정임 그러니까 거의 필요없음 그러니까 안정리할거임

# bcftools reheader --samples names.txt -o /disk2/kb/G400/EB5-IP/EB5-IP.gvcf.gz /disk2/kb/G400/EB5-IP/EB5-IP.gvcf.gz
# tabix -p vcf /disk2/kb/G400/EB5-IP/EB5-IP.gvcf.gz

## import DB 오래걸림
gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport \
    -R ${REFERENCE_FILE} \
    -V ${GVCF} \
    --genomicsdb-workspace-path ${OUT_DIR}/db \
    -L /disk2/kb/sm/T7/total.bed \
    --reader-threads 10 \
    --merge-input-intervals

## Gvcf genotyping
gatk --java-options "-Xmx8g -Xms8g" GenotypeGVCFs \
    -R ${REFERENCE_FILE} \
    -V gendb:/${OUT_DIR}/db\
    -D /knownSite/Homo_sapiens.dbsnp155_newcontig.vcf.gz \
    -O ${OUT_DIR}/hg38.vcf.gz



# bgzip /disk2/kb/G400/vcf/hg38.vcf
# tabix -p vcf /disk2/kb/G400/vcf/hg38.vcf.gz


## SNP과 INDEL을 나눔
gatk --java-options "-Xmx4G" SelectVariants \
            -R /reference/hg38.fa.gz \
            -V %s/hg38.vcf.gz \
            -L /disk2/kb/sm/T7/total.bed \
            --select-type-to-include SNP \
            -O %s/hg38_snps.vcf.gz' % (name, name))

gatk --java-options "-Xmx4G" SelectVariants \
            -R /reference/hg38.fa.gz \
            -V %s/hg38.vcf.gz \
            -L /disk2/kb/sm/T7/total.bed \
            --select-type-to-include INDEL \
            -O %s/hg38_indels.vcf.gz' % (name, name))

## 위치 정보만 가져옴(빠른 일처리를 위함)
gatk --java-options "-Xmx4G" MakeSitesOnlyVcf \
            -I %s/hg38_snps.vcf.gz \
            -O %s/hg38_snps_sitesonly.vcf.gz' % (name, name))

gatk --java-options "-Xmx4G" MakeSitesOnlyVcf \
            -I %s/hg38_indels.vcf.gz \
            -O %s/hg38_indels_sitesonly.vcf.gz' % (name, name))

## 빠르게 필터링함(scoring)
gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
            -R /reference/hg38.fa.gz \
            -V %s/hg38_snps_sitesonly.vcf.gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            --mode SNP \
            --max-gaussians 6 \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /knownSite/hapmap_3.3.hg38.vcf.gz \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 /knownSite/1000G_omni2.5.hg38.vcf.gz \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 /knownSite/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 /knownSite/Homo_sapiens.dbsnp155_newcontig.vcf.gz \
            -O %s/snps.recal \
            --tranches-file %s/snps.tranches' % (name, name, name))

gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
            -R /reference/hg38.fa.gz \
            -V %s/hg38_indels_sitesonly.vcf.gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            --mode INDEL \
            --max-gaussians 4 \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 /knownSite/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 /knownSite/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /knownSite/Homo_sapiens.dbsnp155_newcontig.vcf.gz \
            -O %s/indels.recal \
            --tranches-file %s/indels.tranches' % (name, name, name))

### 적용 후 필터링
gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -R /reference/hg38.fa.gz \
            -V %s/hg38_snps.vcf.gz \
            --recal-file %s/snps.recal \
            --tranches-file %s/snps.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            --mode SNP \
            -O %s/hg38_snps_recalibrated.vcf.gz' % (name, name, name, name))

gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -R /reference/hg38.fa.gz \
            -V %s/hg38_indels.vcf.gz \
            --recal-file %s/indels.recal \
            --tranches-file %s/indels.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            --mode INDEL \
            -O %s/hg38_indels_recalibrated.vcf.gz' % (name, name, name, name))

## 진짜 필터링 
gatk --java-options "-Xmx4G" VariantFiltration \
            -V %s/hg38_snps_recalibrated.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O %s/hg38_snps_filtered.vcf.gz' % (name, name))

gatk --java-options "-Xmx4G" VariantFiltration \
            -V %s/hg38_indels_recalibrated.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O %s/hg38_indels_filtered.vcf.gz' % (name, name))

## snp이랑 indel 합치기
bcftools concat -a %s/hg38_snps_filtered.vcf.gz \
    %s/hg38_indels_filtered.vcf.gz | \
    bgzip > %s/hg38_filtered.vcf.gz' % (name, name, name))

## indexing
tabix -p vcf %s/hg38_filtered.vcf.gz' % (name))

# rm -rf %s/hg38_snps*.vcf.gz* %s/hg38_indels*.vcf.gz* %s/.recal* %s/.tranches %s/_sitesonly.vcf.gz
# rm -rf intervals.list

# done < sample_groups.txt
