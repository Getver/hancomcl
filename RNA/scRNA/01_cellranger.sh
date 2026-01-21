SAMPLE=$1
REF='/usr/local/src/reference/refdata-gex-GRCh38-2024-A/'

cellranger count \
	--id=${SAMPLE} \
	--transcriptome=${REF} \
	--fastqs=/home/data/fastq \
	--sample=${SAMPLE} \
	--create-bam=true

# https://www.embopress.org/doi/full/10.15252/msb.20188746
