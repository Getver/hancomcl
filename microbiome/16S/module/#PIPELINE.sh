PIPELINE () {
    # local INPUT_FILE="$1"
    # local OUTPUT_NAME="$2"
    local LOGFILE="/home/data/logs/${OUTPUT_NAME}_runtime.log"
    local MANIFEST_FILE="$1"
    local OUTPUT_NAME="$2"
    local DOWN_STREAM="${DOWN_STREAM:-false}"

    
    shift 2
    if [[ ${#GROUP_COLUMNS[@]} -eq 0 ]]; then
        GROUP_COLUMNS=("group")
    fi

    # mkdir -p logs
    echo "=== Runtime log for ${OUTPUT_NAME} ===" > "$LOGFILE"

    log_runtime() {
        local step_name="$1"
        shift
        local start_time=$(date +%s)
        if ! "$@"; then
            echo "[ERROR] $step_name 실패" >> "$LOGFILE"
            return 1
        fi
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        printf "[%s] %02d:%02d:%02d\n" \
            "$step_name" $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60)) \
            >> "$LOGFILE"
    }





     if [[ "$SERVICE" == 'Y' ]]; then
        log_runtime "Inspect_metadata" Inspect_metadata
        log_runtime "Denoise_table" Denoise_table
        log_runtime "Taxonomic_analysis" Taxonomic_analysis
        log_runtime "Filtered_noMT" Filtered_noMT
        log_runtime "Calculation_abundance" Calculation_abundance

        echo "[INFO] Y존 서비스 옵션 → blast 단계 실행 중..."
        log_runtime "Ncbi_blast" Ncbi_blast

        python3 /work/script/ncbi_matching.py

        log_runtime "Taxonomic_barplot" Taxonomic_barplot

        python3 /work/script/ASV_check.py
        python3 /work/script/Y/run.py

    elif [[ "$SERVICE" == 'QC' ]]; then
        log_runtime "Inspect_metadata" Inspect_metadata
        log_runtime "Denoise_table" Denoise_table
        log_runtime "Taxonomic_analysis" Taxonomic_analysis
        log_runtime "Filtered_noMT" Filtered_noMT
        log_runtime "Calculation_abundance" Calculation_abundance
        
        log_runtime "Ncbi_blast" Ncbi_blast # ncbi DB 활용하여 blast 해주는 코드

        python3 /work/script/ncbi_matching.py   # ncbi 파일을 유저에 맞게 RESULT_TABLE로 뽑아주는 파일, 수정 필요

        log_runtime "Taxonomic_barplot" Taxonomic_barplot   # ncbi_blast 결과를 이용하여 barplot을 그려주는 파일

        python3 /work/script/ASV_check.py   # qiime2의 silva DB genus 까지의 분석 결과와 ncbi의 16S_rna 분석 결과를 비교하여 최상의 tsv를 뽑아주는 파일

        log_runtime "Make_tree" Make_tree   # abundance, taxonomy 파일을 참조하여 phylogenic 분석하여 tree 분석해주는 코드

        log_runtime "Diversity_AB" Diversity_AB # alpha, beta diversity 분석해주는 코드

        python3 /work/script/transfer_krona.py  # krona input file 만들어주는 코드
        ktImportText taxonomy/krona_input.tsv -o taxonomy/krona.html    # 잘 돌아가는데 그림 안보임, html 문제같음, 확인 해야됨

        Rscript /work/script/plot_tree.R    #  phylogenic tree 그려주는 파일, 코드 수정해야할듯

        if [[ "$DOWN_STREAM" == true && ${#GROUP_COLUMNS[@]} -gt 0 ]]; then
            echo "[INFO] Analysis_group columns: ${GROUP_COLUMNS[*]}"
            log_runtime "Analysis_group" Analysis_group "${GROUP_COLUMNS[@]}"
            Rscript /work/script/ANCOM.R "${GROUP_COLUMNS[@]}"
        fi

    else
        # SERVICE 값이 비어있거나 다른 값일 때 기본 QC 파이프라인
        echo "[INFO] SERVICE 값이 비어 있거나 알 수 없는 값 → 기본 QC 파이프라인 실행"
        log_runtime "Inspect_metadata" Inspect_metadata
        log_runtime "Denoise_table" Denoise_table
        log_runtime "Taxonomic_analysis" Taxonomic_analysis
        log_runtime "Filtered_noMT" Filtered_noMT
        log_runtime "Calculation_abundance" Calculation_abundance
        
        log_runtime "Ncbi_blast" Ncbi_blast

        python3 /work/script/ncbi_matching.py

        log_runtime "Taxonomic_barplot" Taxonomic_barplot
        log_runtime "Make_tree" Make_tree

        log_runtime "Diversity_AB" Diversity_AB

        python3 /work/script/transfer_krona.py
        ktImportText taxonomy/krona_input.tsv -o taxonomy/krona.html

        Rscript /work/script/plot_tree.R

        if [[ "$DOWN_STREAM" == true && ${#GROUP_COLUMNS[@]} -gt 0 ]]; then
            echo "[INFO] Analysis_group columns: ${GROUP_COLUMNS[*]}"
            log_runtime "Analysis_group" Analysis_group "${GROUP_COLUMNS[@]}"
            Rscript /work/script/ANCOM.R "${GROUP_COLUMNS[@]}"
        fi
    fi


    echo "=== Pipeline finished for ${OUTPUT_NAME} ===" >> "$LOGFILE"

}



