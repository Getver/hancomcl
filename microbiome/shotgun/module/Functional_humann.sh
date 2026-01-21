## Functional profiling : humann
Functional_humann () {
    local SAMPLE="$1"

    mkdir -p humann

    humann \
        --input host_removed/${SAMPLE}_merged.fastq.gz \
        --output humann/${SAMPLE} \
        --nucleotide-database /work/DB/humann_db/chocophlan \
        --protein-database /work/DB/humann_db/uniref \
        --threads 70 \
        --remove-temp-output \
        --taxonomic-profile metaphlan/${SAMPLE}_profile.txt
    
    # CPM으로 정규화
    humann_renorm_table \
        --input humann/${SAMPLE}/${SAMPLE}_merged_genefamilies.tsv \
        --output humann/${SAMPLE}/${SAMPLE}_merged_genefamilies_cpm.tsv \
        --units cpm \
        --update-snames

    # KO 단위로 변환
    humann_regroup_table \
        --input humann/${SAMPLE}/${SAMPLE}_merged_genefamilies.tsv \
        --custom /work/DB/humann_db/utility_mapping/map_ko_uniref90.txt.gz \
        --output humann/${SAMPLE}/${SAMPLE}_merged_genefamilies_ko.tsv
}
# real    71m42.390s
# user    751m56.097s
# sys     8m11.473s