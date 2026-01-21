
## GPU
# https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_germline.html#man-germline

INPUT_FASTQ_1='T7/ztron/ztron_project_id/R2100610230006/E200014516/L01/E200014516_L01/23110600120S10B1A100000A00_R1.fastq.gz'
INPUT_FASTQ_2='T7/ztron/ztron_project_id/R2100610230006/E200014516/L01/E200014516_L01/23110600120S10B1A100000A00_R2.fastq.gz'
OUTPUT_BAM='sm/MEI/T7test/vcf2/23110600120S10B1A100000A00.bam'
OUTPUT_VCF='sm/MEI/T7test/vcf2/23110600120S10B1A100000A00.vcf'
OUT_RECAL_FILE='sm/MEI/T7test/vcf2/23110600120S10B1A100000A00.recal'
REFERENCE_FILE='sm/WGS/reference/hg38.fa.gz'
KNOWN_SITES_FILE='references/hg38/UCSC/DNA/bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'


INPUT_DIR='/disk2'
OUTPUT_DIR='/disk2'


# This command assumes all the inputs are in INPUT_DIR and all the outputs go to OUTPUT_DIR.
docker run --rm --gpus all --volume ${INPUT_DIR}:/workdir --volume ${OUTPUT_DIR}:/outputdir \
    --workdir /workdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
    pbrun germline \
    --ref /workdir/${REFERENCE_FILE} \
    --in-fq /workdir/${INPUT_FASTQ_1} /workdir/${INPUT_FASTQ_2} \
    --knownSites /workdir/${KNOWN_SITES_FILE} \
    --out-bam /outputdir/${OUTPUT_BAM} \
    --out-variants /outputdir/${OUTPUT_VCF} \
    --out-recal-file /outputdir/${OUT_RECAL_FILE} \
    --tmp-dir /workdir/sm/MEI/T7test/tmpdir

