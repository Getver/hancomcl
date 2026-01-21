Result_merge () {
    ## metaphlan result merge
    ## 샘플 taxonomy annotation 분석 결과 합치기
    /usr/local/bin/merge_metaphlan_tables.py metaphlan/*_profile.txt > metaphlan/merged_abundance_table.txt
    ## SGB 이름 GTDB로 변환(metaphlan4 부터 도입된 기능인데.. 나는 humann과의 호환때메 3밖에 못써서 git에서 퍼옴...... 인데 형식 달라져서 에러남... 걍 못함... ㅠ)
    # /work/module/sgb_to_gtdb/sgb_to_gtdb_profile.py -i metaphlan/merged_abundance_table.txt -o metaphlan/merged_abundance_table_gtdb.txt

    ## humann result merge
    humann_join_tables -i humann -o humann/merged_subset_genefamilies.tsv --file_name genefamilies
}