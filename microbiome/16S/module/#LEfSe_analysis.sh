LEfSe_analysis() {
    local GROUP_COLUMNS=("$@")
    [[ ${#GROUP_COLUMNS[@]} -eq 0 ]] && GROUP_COLUMNS=("group")
    # GROUP_COLUMNS=("group" "day")

    mkdir -p LEfSe

    python3 /work/script/make_LEfSe.py

    for GROUP_COLUMN in "${GROUP_COLUMNS[@]}"; do
        # 컬럼별 폴더 생성
        local COL_DIR="LEfSe/$GROUP_COLUMN"
        mkdir -p "$COL_DIR"

        local LEFSE_INPUT="LEfSe/lefse_input.in"
        local LEFSE_RESULTS="LEfSe/lefse_results.res"
        local LEFSE_PLOT="LEfSe/lefse_plot.png"

        # LEfSe GitHub Python 스크립트 경로
        local LEFSE_DIR="/work/tool/lefse"

        # 포맷 변환
        TOTAL_COL=$(python3 -c "import pandas as pd; import sys; df = pd.read_csv('LEfSe/lefse_input_combined.tsv', sep='\t'); print(df.shape[0]-1)")
        python3 /work/tool/lefse/lefse_format_input.py LEfSe/lefse_input_combined.tsv LEfSe/lefse_input.in -c 2 -u 1 -o "$TOTAL_COL"

        # LEfSe 실행
        python3 /work/tool/lefse/lefse_run.py "$LEFSE_INPUT" "$LEFSE_RESULTS"

        # 플롯 생성
        python3 /work/tool/lefse/lefse_plot_res.py "$LEFSE_RESULTS" "$LEFSE_PLOT" --format png

        # echo "[INFO] LEfSe 분석 완료 → $COL_DIR"
    done
}
