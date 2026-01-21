#!/bin/bash
# microbiome shotgun sequencing read-based annotation PIPELINE
# Inside Docker
source /work/module/Sanitize.sh
source /work/module/PIPELINE.sh
source /work/module/Inspect_metadata.sh
source /work/module/Denoise_table.sh
source /work/module/Taxonomic_analysis.sh
source /work/module/Filtered_noMT.sh
source /work/module/Taxonomic_barplot.sh
source /work/module/Make_tree.sh
source /work/module/Diversity_AB.sh
source /work/module/Calculation_abundance.sh
source /work/module/Ncbi_blast.sh
source /work/module/Analysis_group.sh

set -euo pipefail


# ---------------------------
# 0) 환경 변수 (Docker 내부)
# ---------------------------
export WORK_DIR="/home/data"
export RAW_DATA="/work/raw_data"
export MANIFEST="/work/input_data.tsv"      # 항상 여기 생성
export METADATA="/work/meta_data.tsv"       # down-stream 시 생성
export DOWN_STREAM="${DOWN_STREAM:-false}"
export SERVICE="${SERVICE:-}"
export META_DATA="${META_DATA:-}"           # Docker 내부 경로
export DOWN_STREAM_COLUMNS="${DOWN_STREAM_COLUMNS:-}"





OUTPUT_NAME="$1"

echo "[INFO] OUTPUT_NAME: $OUTPUT_NAME"
[[ -n "$META_DATA" ]] && echo "[INFO] META_DATA (docker 내부): $META_DATA"
echo "[INFO] SERVICE: ${SERVICE:-"(none)"}"
echo "[INFO] DOWN_STREAM: $DOWN_STREAM"
[[ -n "$DOWN_STREAM_COLUMNS" ]] && echo "[INFO] DOWN_STREAM_COLUMNS: $DOWN_STREAM_COLUMNS"

# ---------------------------
# 1) sample-id 추출
# ---------------------------
if [[ -n "$META_DATA" ]]; then
    echo "[INFO] META_DATA 제공됨 → sample-id 추출: $META_DATA"

    if [[ ! -f "$META_DATA" ]]; then
        echo "[ERROR] META_DATA 파일을 찾을 수 없음: $META_DATA"
        exit 1
    fi

    SAMPLE_IDS=$(tail -n +2 "$META_DATA" | cut -f1 | xargs -n1 basename)
else
    echo "[INFO] META_DATA 없음 → RAW_DATA 기반 sample-id 자동 생성"

    FILE_LIST=$(ls "$RAW_DATA"/*_R1*fastq* "$RAW_DATA"/*_1.fastq* "$RAW_DATA"/*.R1.* 2>/dev/null || true)
    if [[ -z "$FILE_LIST" ]]; then
        echo "[ERROR] /work/raw_data 에서 R1 FASTQ 파일을 찾을 수 없음"
        exit 1
    fi

    SAMPLE_IDS=""
    for f in $FILE_LIST; do
        base=$(basename "$f")
        case "$base" in
            *_R1_*) sample="${base%%_R1_*}" ;;
            *_R1.*) sample="${base%%_R1.*}" ;;
            *.R1.*) sample="${base%%.R1.*}" ;;
            *_1.fastq*) sample="${base%%_1.fastq*}" ;;
            *) echo "[WARN] FASTQ 이름 패턴 불명 → 스킵: $base"; continue ;;
        esac
        SAMPLE_IDS="$SAMPLE_IDS $sample"
    done

    SAMPLE_IDS=$(echo "$SAMPLE_IDS" | tr ' ' '\n' | sort -u)
fi

# ---------------------------
# 2) DOWN_STREAM → meta_data.tsv 생성
# ---------------------------
if [[ "$DOWN_STREAM" == true ]]; then
    if [[ -z "$META_DATA" ]]; then
        echo "[ERROR] --down-stream 사용 시 --meta-data 옵션 필수"
        exit 1
    fi
    echo "[INFO] --down-stream 활성화 → meta_data.tsv 생성 중..."
    sanitize_metadata "$META_DATA" > "$METADATA"
    echo "[INFO] meta_data.tsv 생성 완료 → $(wc -l < "$METADATA") lines"
fi

if [[ -z "$SAMPLE_IDS" ]]; then
    echo "[ERROR] sample-id 추출 실패"
    exit 1
fi

# ---------------------------
# 3) MANIFEST 생성
# ---------------------------
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

for SAMPLE in $SAMPLE_IDS; do
    R1=$(ls "$RAW_DATA"/"$SAMPLE"*R1* "$RAW_DATA"/"$SAMPLE"*1.fastq* "$RAW_DATA"/"$SAMPLE"*.R1.* 2>/dev/null | head -n1 || true)
    R2=$(ls "$RAW_DATA"/"$SAMPLE"*R2* "$RAW_DATA"/"$SAMPLE"*2.fastq* "$RAW_DATA"/"$SAMPLE"*.R2.* 2>/dev/null | head -n1 || true)

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo -e "${SAMPLE}\t${R1}\t${R2}" >> "$MANIFEST"
    else
        echo "[WARN] R1 또는 R2 FASTQ 없음 → 스킵: $SAMPLE"
    fi
done

echo "[INFO] MANIFEST 생성 완료 → $MANIFEST"

# ---------------------------
# 4) Downstream 컬럼 배열 변환
# ---------------------------
if [[ -n "$DOWN_STREAM_COLUMNS" ]]; then
    IFS=',' read -ra GROUP_COLUMNS <<< "$DOWN_STREAM_COLUMNS"
else
    GROUP_COLUMNS=("group")  # 기본값
fi
echo "[INFO] GROUP_COLUMNS=${GROUP_COLUMNS[*]}"

# # ---------------------------
# # 4-1) Sequencing Platform 판별 (Illumina vs Nanopore)
# # ---------------------------
# echo "[INFO] Sequencing Platform 자동 판별 중..."

# # MANIFEST의 첫 번째 데이터 라인에서 Forward 파일 경로 추출
# FIRST_FASTQ=$(awk 'NR==2 {print $2}' "$MANIFEST")

# if [[ -z "$FIRST_FASTQ" || ! -f "$FIRST_FASTQ" ]]; then
#     echo "[ERROR] 판별할 FASTQ 파일을 찾을 수 없습니다: $FIRST_FASTQ"
#     exit 1
# fi

# # 압축 파일(*.gz) 대응
# if [[ "$FIRST_FASTQ" == *.gz ]]; then
#     PEEK_CMD="gzip -cd"
# else
#     PEEK_CMD="cat"
# fi

# # 처음 100개 리드(400라인)의 서열 길이 평균 계산
# # awk: NR%4==2 (서열 라인)일 때 길이 합산, 400라인 넘으면 종료
# # pipefail 에러 방지를 위해 '|| true' 추가
# AVG_READ_LEN=$($PEEK_CMD "$FIRST_FASTQ" | head -n 400 | awk '{if(NR%4==2) {sum+=length($0); cnt++}} END {if (cnt>0) print int(sum/cnt); else print 0}' || true)

# # 판별 로직 (Threshold: 600bp)
# # - Illumina: 보통 150bp (PE) ~ 300bp (MiSeq) -> 600bp 미만
# # - Nanopore: 보통 1kb 이상 -> 600bp 이상
# if [[ "$AVG_READ_LEN" -gt 600 ]]; then
#     export PLATFORM="Nanopore"
# else
#     export PLATFORM="Illumina"
# fi

# echo "[INFO] 첫 100개 리드 평균 길이: ${AVG_READ_LEN} bp"
# echo "[INFO] 감지된 플랫폼: $PLATFORM"

# ---------------------------
# 5) 파이프라인 실행
# ---------------------------
PIPELINE "$MANIFEST" "$OUTPUT_NAME" "${GROUP_COLUMNS[@]}"c
