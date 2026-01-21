# import schedule
import time
import glob, os

# docker run --rm -it -v /disk0/kb/kobic/vep:/cache -v /disk0/kb/vep:/plugins vc:latest

# vep --cache --dir_cache /cache --dir_plugins /plugins --offline -i <input> -o <output>

# /disk0/sm/WES/tool/SNPlocs.Hsapiens.dbSNP151.GRCh38/inst/extdata/rowids.rds

def run_bwa(threads, mech_name, fc_num, lane_num, ref, out_dir, knownSite_dir, adaptor_num):
    for AD in adaptor_num :
        path = '/disk2/' + '/'.join([mech_name, fc_num, lane_num])
        fName = '_'.join([fc_num, lane_num, AD])
        fList = glob.glob(path + '/' + fName + '_*.gz')
        #
        ## fastqc
        os.system('fastqc -o %s %s %s' %(out_dir+'/fastqc', fList[0], fList[1]))
        #
        ## read group 지정(꼭 저렇게 \\t 이렇게 써야됨, ' 있어야함)
        rg = "'@RG\\tID:%s\\tSM:%s\\tPL:DNBSEQ\\tPU:flowcell_barcode\'" %(AD, adaptor_num[AD])
        name = '%s%s/%s' % (out_dir, adaptor_num[AD], adaptor_num[AD])
        #
        ## bam 만들기
        os.system('bwa mem -t %s %s %s %s -R %s | samtools view -h -b -o %s_hg38.bam' %(threads, ref, fList[0], fList[1], rg, name))
        #
        ## sortsam 써서 sort_bam 만들기
        os.system('gatk --java-options "-Xmx4G" SortSam \
                        -I %s_hg38.bam \
                        -O %s_hg38_sorted.bam -SO coordinate \
                        --VALIDATION_STRINGENCY SILENT \
                        --MAX_RECORDS_IN_RAM 5000000' % (name, name))
        os.system('samtools index %s_hg38_sorted.bam' % (name))
        #
        ## on-target 얻고싶을때 쓰기
        os.system('java -jar /disk2/sm/cancer/tool/picard.jar CollectHsMetrics \
                        --INPUT %s_hg38_sorted.bam \
                        --OUTPUT %s_hsmetrics.txt \
                        --BAIT_INTERVALS /disk2/kb/G400/target.interval_list \
                        --TARGET_INTERVALS /disk2/kb/G400/target.interval_list' % (name, name))
        #
        ## qualimap(sorted 하고나서 해야함)
        os.system('qualimap --java-mem-size=80G bamqc \
                            -bam %s_hg38_sorted.bam \
                            --outfile %s' %(name, adaptor_num[AD]))
        #
        ## markduplicates
        os.system('gatk --java-options "-Xmx4G" MarkDuplicates \
                        -I %s_hg38_sorted.bam \
                        -O %s_hg38_dup.bam \
                        -M %s_hg38_dup.metrics \
                        --VALIDATION_STRINGENCY SILENT \
                        --MAX_RECORDS_IN_RAM 5000000' % (name, name, name))
        os.system('samtools index %s_hg38_dup.bam')
        #
        ## BaseRecalibrator 다음부터는 하지 않을 스텝임(너무 오래걸림)
        os.system('gatk --java-options "-Xmx4G" BaseRecalibrator -R %s -I %s_hg38_dup.bam \
                        --known-sites %s/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                        --known-sites %s/hapmap_3.3.hg38.vcf.gz \
                        --known-sites %s/1000G_omni2.5.hg38.vcf.gz \
                        --known-sites %s/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
                        --known-sites %s/Homo_sapiens_assembly38.known_indels.vcf.gz \
                        --known-sites %s/Homo_sapiens_hg38.dbsnp155_final.vcf.gz \
                        -O %s_hg38_dup.recal' % (ref, name, knownSite_dir, knownSite_dir, knownSite_dir, knownSite_dir, knownSite_dir, knownSite_dir, name))
        os.system('gatk --java-options "-Xmx4G" ApplyBQSR -R %s -I %s_hg38_dup.bam --bqsr-recal-file %s_hg38_dup.recal -O %s_hg38_recal.bam \
                        --create-output-bam-index false' % (ref, name, name, name))
        os.system('samtools index %s_hg38_recal.bam' % (name))
        #
        ## multiplemetrics 얻고싶을때 쓰기
        os.system('gatk --java-options "-Xmx4G" CollectMultipleMetrics \
                        -R %s \
                        -I %s_hg38_recal.bam \
                        -O %s_hg38_aligned \
                        -AS true \
                        --PROGRAM "CollectAlignmentSummaryMetrics" \
                        --PROGRAM "CollectInsertSizeMetrics" \
                        --PROGRAM "CollectSequencingArtifactMetrics" \
                        --PROGRAM "CollectGcBiasMetrics" \
                        --PROGRAM "QualityScoreDistribution"' % (ref, name, name))
        #
        ## haplotypecaller
        os.system('gatk --java-options "-Xms8G" HaplotypeCaller \
                        -R %s \
                        -I %s_hg38_recal.bam \
                        -O %s_hg38.gvcf.gz \
                        -D %s/Homo_sapiens_hg38.dbsnp155_final.vcf.gz \
                        -ERC GVCF' % (ref, name, name, knownSite_dir))


def run_gatk(ref, knownSite_dir, name):
    ## import DB
    os.system('gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport \
                    -R %s \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/EB4-IP_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/EB4-2W_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/EB4-4W_hg38.gvcf.gz \
                    -V /disk2/kb/G400/EB5-IP/EB5-IP.gvcf.gz \
                    -V /disk2/kb/G400/EB5-2W/EB5-2W.gvcf.gz \
                    -V /disk2/kb/G400/EB5-4W/EB5-4W.gvcf.gz \
                    -V /disk2/kb/G400/EB9-IP/EB9-IP.gvcf.gz \
                    -V /disk2/kb/G400/EB9-2W/EB9-2W.gvcf.gz \
                    -V /disk2/kb/G400/EB9-4W/EB9-4W.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB2-IP_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB2-2W_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB2-4W_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB3-IP_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB3-2W_hg38.gvcf.gz \
                    -V /disk2/kb/exome_sequencing/all_samples_analysis/gvcfs/LB3-4W_hg38.gvcf.gz \
                    -V /disk2/kb/G400/LB4-IP/LB4-IP.gvcf.gz \
                    -V /disk2/kb/G400/LB4-2W/LB4-2W.gvcf.gz \
                    -V /disk2/kb/G400/LB4-4W/LB4-4W.gvcf.gz \
                    --genomicsdb-workspace-path /disk2/kb/G400/db \
                    -L /disk2/kb/G400/V5_V8_intersect_hg38.bed \
                    --reader-threads 9 \
                    --merge-input-intervals' % (ref))
    #
    ## Gvcf genotyping
    os.system('gatk --java-options "-Xmx8g -Xms8g" GenotypeGVCFs \
                    -R %s \
                    -V gendb://db\
                    -D %s/Homo_sapiens_hg38.dbsnp155_newcontig.vcf.gz \
                    -O %s/G400_hg38.vcf.gz' % (ref, knownSite_dir, name))
    #
    ## SNP과 INDEL을 나눔
    os.system('gatk --java-options "-Xmx4G" SelectVariants \
                    -R %s \
                    -V %s/G400_hg38.vcf.gz \
                    -L /disk2/kb/G400/V5_V8_intersect_hg38.bed \
                    --select-type-to-include SNP \
                    -O %s/G400_hg38_snps.vcf.gz' % (ref, name, name))
    #
    os.system('gatk --java-options "-Xmx4G" SelectVariants \
                    -R %s \
                    -V %s/G400_hg38.vcf.gz \
                    -L /disk2/kb/G400/V5_V8_intersect_hg38.bed \
                    --select-type-to-include INDEL \
                    -O %s/G400_hg38_indels.vcf.gz' % (ref, name, name))
    #
    ## 위치 정보만 가져옴(빠른 일처리를 위함)
    os.system('gatk --java-options "-Xmx4G" MakeSitesOnlyVcf \
                    -I %s/G400_hg38_snps.vcf.gz \
                    -O %s/G400_hg38_snps_sitesonly.vcf.gz' % (name, name))

    os.system('gatk --java-options "-Xmx4G" MakeSitesOnlyVcf \
                    -I %s/G400_hg38_indels.vcf.gz \
                    -O %s/G400_hg38_indels_sitesonly.vcf.gz' % (name, name))
    #
    ## 빠르게 필터링함(scoring)
    os.system('gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
                    -R %s \
                    -V %s/G400_hg38_snps_sitesonly.vcf.gz \
                    --trust-all-polymorphic \
                    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
                    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
                    --mode SNP \
                    --max-gaussians 6 \
                    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s/hapmap_3.3.hg38.vcf.gz \
                    --resource:omni,known=false,training=true,truth=true,prior=12.0 %s/1000G_omni2.5.hg38.vcf.gz \
                    --resource:1000G,known=false,training=true,truth=false,prior=10.0 %s/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
                    --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 %s/Homo_sapiens_hg38.dbsnp155_newcontig.vcf.gz \
                    -O %s/G400_snps.recal \
                    --tranches-file %s/G400_snps.tranches' % (ref, name, knownSite_dir, knownSite_dir, knownSite_dir, knownSite_dir, name, name))
    #
    os.system('gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
                    -R %s \
                    -V %s/G400_hg38_indels_sitesonly.vcf.gz \
                    --trust-all-polymorphic \
                    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
                    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
                    --mode INDEL \
                    --max-gaussians 4 \
                    --resource:mills,known=false,training=true,truth=true,prior=12.0 %s/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                    --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 %s/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
                    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s/Homo_sapiens_hg38.dbsnp155_newcontig.vcf.gz \
                    -O %s/G400_indels.recal \
                    --tranches-file %s/G400_indels.tranches' % (ref, name, knownSite_dir, knownSite_dir, knownSite_dir, name, name))
    #
    ### 적용 후 필터링
    os.system('gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
                    -R %s \
                    -V %s/G400_hg38_snps.vcf.gz \
                    --recal-file %s/G400_snps.recal \
                    --tranches-file %s/G400_snps.tranches \
                    --truth-sensitivity-filter-level 99.7 \
                    --create-output-variant-index true \
                    --mode SNP \
                    -O %s/G400_hg38_snps_recalibrated.vcf.gz' % (ref, name, name, name, name))
    #
    os.system('gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
                    -R %s \
                    -V %s/G400_hg38_indels.vcf.gz \
                    --recal-file %s/G400_indels.recal \
                    --tranches-file %s/G400_indels.tranches \
                    --truth-sensitivity-filter-level 99.7 \
                    --create-output-variant-index true \
                    --mode INDEL \
                    -O %s/G400_hg38_indels_recalibrated.vcf.gz' % (ref, name, name, name, name))
    #
    ## 진짜 필터링 
    os.system('gatk --java-options "-Xmx4G" VariantFiltration \
                    -V %s/G400_hg38_snps_recalibrated.vcf.gz \
                    -filter "QD < 2.0" --filter-name "QD2" \
                    -filter "QUAL < 30.0" --filter-name "QUAL30" \
                    -filter "SOR > 3.0" --filter-name "SOR3" \
                    -filter "FS > 60.0" --filter-name "FS60" \
                    -filter "MQ < 40.0" --filter-name "MQ40" \
                    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                    -O %s/G400_hg38_snps_filtered.vcf.gz' % (name, name))
    #
    os.system('gatk --java-options "-Xmx4G" VariantFiltration \
                    -V %s/G400_hg38_indels_recalibrated.vcf.gz \
                    -filter "QD < 2.0" --filter-name "QD2" \
                    -filter "QUAL < 30.0" --filter-name "QUAL30" \
                    -filter "FS > 200.0" --filter-name "FS200" \
                    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
                    -O %s/G400_hg38_indels_filtered.vcf.gz' % (name, name))
    #
    ## snp이랑 indel 합치기
    os.system('bcftools concat -a %s/G400_hg38_snps_filtered.vcf.gz \
            %s/G400_hg38_indels_filtered.vcf.gz | \
            bgzip > %s/G400_hg38_filtered.vcf.gz' % (name, name, name))
    #
    ## indexing
    os.system('tabix -p vcf %s/G400_hg38_filtered.vcf.gz')


def run_calling():
    threads = '32'
    mech_name = 'C2130410230003'
    fc_num = 'F350023255'
    lane_num = 'L01'
    ref = '/disk2/references/hg38/UCSC/DNA/all_chromosomes/hg38.fa.gz'
    out_dir = '/disk2/kb/G400/'
    knownSite_dir="/disk2/kb/exome_sequencing"
    ## G400전용 adaptor 넣기
    adaptor_num = {'1':'EB5-IP', 
                   '2':'EB5-2W', 
                   '3':'EB5-4W', 
                   '4':'EB9-IP', 
                   '13':'EB9-2W', 
                   '14':'EB9-4W', 
                   '15':'LB4-IP',
                   '16':'LB4-2W', 
                   '97':'LB4-4W', 
                   '98':'CLP00046', 
                   '99':'CLP00050', 
                   '100':'CLP00051', 
                   '101':'CLP00052', 
                   '102':'CLP00057', 
                   '103':'CLP00052_2', 
                   '104':'CLP00057_2'}
    #
    gatk_name = out_dir+'vcf'
    run_bwa(threads, mech_name, fc_num, lane_num, ref, out_dir, knownSite_dir, adaptor_num)
    run_gatk(ref, knownSite_dir, gatk_name)
    

run_calling()


# def sumin():
    # print('sumin')

# schedule.every().tuesday.at('17:20').do(sumin)
# schedule.every(5).seconds.do(sumin)
# schedule.every().thursday.at('6:00').do(run_calling)
# schedule.clear(tag=None)
# while True:
#     schedule.run_pending()
#     time.sleep(1)