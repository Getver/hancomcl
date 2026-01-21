# https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#emptydrop-like-filtering
# https://docs.tinybio.cloud/docs/starsolo-tutoral
# https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

# SAMPLE=$1
# SAMPLE='pbmc8k_S1'
SAMPLE='sc5p_v2_hs_PBMC_10k_5gex'
REF='/usr/local/src/reference/refdata-gex-GRCh38-2024-A/star'
# --readFilesIn fastq/pbmc8k_S1_L007_R1_001.fastq.gz,fastq/pbmc8k_S1_L008_R1_001.fastq.gz fastq/pbmc8k_S1_L007_R2_001.fastq.gz,fastq/pbmc8k_S1_L008_R2_001.fastq.gz \

mkdir output

/usr/local/src/STAR-2.7.1a/bin/Linux_x86_64/STAR \
     --runThreadN 10 \
     --genomeDir ${REF} \
     --readFilesIn fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R1_001.fastq.gz,fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L002_R1_001.fastq.gz fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R2_001.fastq.gz,fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L002_R2_001.fastq.gz \
     --soloType Droplet \
     --soloCBwhitelist None \
     --soloCBlen 16 \
     --soloBarcodeReadLength 98 \
     --readFilesCommand zcat \
     --outFileNamePrefix output/${SAMPLE} \
     --outSAMtype BAM SortedByCoordinate \
     --soloUMIlen 12 \
     --limitOutSJcollapsed 3000000 \
     --limitOutSJoneRead 10000 \
     --limitIObufferSize 100000000000 \
     --genomeChrBinNbits 11






SAMPLE='pbmc8k_S1'
REF='/usr/local/src/reference/refdata-gex-GRCh38-2024-A/star'

--readFilesIn fastq/pbmc8k_S1_L007_R2_001.fastq.gz,fastq/pbmc8k_S1_L008_R2_001.fastq.gz fastq/pbmc8k_S1_L007_R1_001.fastq.gz,fastq/pbmc8k_S1_L008_R1_001.fastq.gz \




SAMPLE='sc5p_v2_hs_PBMC_10k_5gex'
REF='/usr/local/src/reference/refdata-gex-GRCh38-2024-A/star'

mkdir STARsolo33
/disk0/sm/TEST/TEST_scRNA/tool/STAR-2.7.1a/bin/Linux_x86_64/STAR \
     --runThreadN 10 \
     --genomeDir /disk0/sm/TEST/TEST_scRNA/docker/reference/refdata-gex-GRCh38-2024-A/star \
     --readFilesIn fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R2_001.fastq.gz,fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L002_R2_001.fastq.gz fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R1_001.fastq.gz,fastq/sc5p_v2_hs_PBMC_10k_5gex_S1_L002_R1_001.fastq.gz \
     --soloType Droplet \
     --soloCBwhitelist None \
     --soloCBlen 16 \
     --soloUMIlen 10 \
     --soloBarcodeReadLength 26 \
     --readFilesCommand zcat \
     --outFileNamePrefix STARsolo33/${SAMPLE} \
     --outSAMtype BAM SortedByCoordinate \
     --limitOutSJcollapsed 3000000 \
     --limitOutSJoneRead 10000 \
     --outSAMmultNmax 1 
     \
     --outSAMprimaryFlag AllBestScore \

22 : 
     --outSAMmultNmax 1 
     --outSAMprimaryFlag AllBestScore \
33 : 
     --outSAMprimaryFlag AllBestScore \


# --soloType : starsolo 알고리즘이 켜지도록 함 \
# --soloCBwhitelist : 있는게 좋음 https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist
# --soloCBlen  \
# --soloBarcodeReadLength  \
# --readFilesCommand  \
# --outFileNamePrefix 
# --soloUMIlen : 기본 바코드 길이
# --genomeChrBinNbits 
# --outSAMprimaryFlag 
# --outSAMmultNmax
# --limitOutSJcollapsed 
# --limitOutSJoneRead 

# --winAnchorMultimapNmax 
# --outFilterMultimapNmax

# Feb 04 13:32:24 ..... started STAR run
# Feb 04 13:32:24 ..... loading genome
# Feb 04 13:32:31 ..... started mapping
# Feb 04 16:17:52 ..... finished mapping
# Feb 04 16:17:54 ..... started sorting BAM
# Feb 04 16:27:10 ..... finished successfully


