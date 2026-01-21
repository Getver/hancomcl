SAMPLE=$1

bash /usr/local/src/code/01_cellranger.sh ${SAMPLE}
Rscript /usr/local/src/code/02_seurat.R ${SAMPLE}

sudo chmod -R 777 /home/data/${SAMPLE}
# bash /usr/local/src/code/00_run.sh sc5p_v2_hs_PBMC_10k_5gex
