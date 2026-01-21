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
export WORKDIR="/home/data"
export RAW_DATA="/work/raw_data"
export MANIFEST="/work/input_data.tsv"    # 항상 이 이름
export METADATA="/work/meta_data.tsv"     # 항상 이 이름
export DOWN_STREAM="${DOWN_STREAM:-false}"
export SERVICE="${SERVICE:-}"
export RAW_METADATA="${RAW_METADATA:-}"
export DOWN_STREAM_COLUMNS="${DOWN_STREAM_COLUMNS:-}"


OUTPUT_NAME="$1"

MANIFEST="$WORKDIR/manifest.tsv"
# METADATA="$WORKDIR/meta_data.tsv"

TMP_MANIFEST_PE="$WORKDIR/.manifest_pe.tmp"
TMP_MANIFEST_SE="$WORKDIR/.manifest_se.tmp"

echo "[INFO] RAW_DATA=$RAW_DATA"
echo "[INFO] OUTPUT_NAME=$OUTPUT_NAME"

: "${RAW_DATA:?RAW_DATA not set}"
: "${WORKDIR:?WORKDIR not set}"

# ---------------------------
# 1) sample-id 수집
# ---------------------------
if [[ -n "$RAW_METADATA" ]]; then
    echo "[INFO] META_DATA 기반 sample-id 추출"
    SAMPLE_IDS=$(tail -n +2 "$RAW_METADATA" | cut -f1 | xargs -n1 basename)
else
    echo "[INFO] FASTQ 파일명 기반 sample-id 추출"
    SAMPLE_IDS=$(ls "$RAW_DATA"/*.fastq* | xargs -n1 basename | sed 's/\..*//' | sort -u)
fi

[[ -z "$SAMPLE_IDS" ]] && { echo "[ERROR] sample-id 없음"; exit 1; }

# temp 초기화
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$TMP_MANIFEST_PE"
echo -e "sample-id\tabsolute-filepath" > "$TMP_MANIFEST_SE"

# ---------------------------
# 2) FASTQ 판별
# ---------------------------
for SAMPLE in $SAMPLE_IDS; do
    R1=$(ls "$RAW_DATA"/"$SAMPLE"*R1* "$RAW_DATA"/"$SAMPLE"*1.fastq* 2>/dev/null | head -n1 || true)
    R2=$(ls "$RAW_DATA"/"$SAMPLE"*R2* "$RAW_DATA"/"$SAMPLE"*2.fastq* 2>/dev/null | head -n1 || true)

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo -e "$SAMPLE\t$R1\t$R2" >> "$TMP_MANIFEST_PE"
    else
        SE_FASTQ=$(ls "$RAW_DATA"/"$SAMPLE".fastq* 2>/dev/null | head -n1 || true)

        if [[ -f "$SE_FASTQ" ]]; then
            echo -e "$SAMPLE\t$SE_FASTQ" >> "$TMP_MANIFEST_SE"
        fi
    fi
done

# ---------------------------
# 3) paired / single 결정
# ---------------------------
PE_COUNT=$(($(wc -l < "$TMP_MANIFEST_PE") - 1))
SE_COUNT=$(($(wc -l < "$TMP_MANIFEST_SE") - 1))

echo "[INFO] PE_COUNT=$PE_COUNT SE_COUNT=$SE_COUNT"

if [[ "$PE_COUNT" -gt 0 && "$SE_COUNT" -gt 0 ]]; then
    echo "[ERROR] paired + single 혼합 데이터"
    exit 1
elif [[ "$PE_COUNT" -gt 0 ]]; then
    cp "$TMP_MANIFEST_PE" "$MANIFEST"
    SEQ_MODE="paired"
elif [[ "$SE_COUNT" -gt 0 ]]; then
    cp "$TMP_MANIFEST_SE" "$MANIFEST"
    SEQ_MODE="single"
else
    echo "[ERROR] FASTQ 인식 실패"
    exit 1
fi

echo "[INFO] SEQ_MODE=$SEQ_MODE"
echo "[INFO] MANIFEST=$MANIFEST"

# ---------------------------
# 4) meta_data.tsv
# ---------------------------
if [[ "$DOWN_STREAM" == true ]]; then
    [[ -z "$RAW_METADATA" ]] && { echo "[ERROR] --down-stream 시 META_DATA 필수"; exit 1; }
    sanitize_metadata "$RAW_METADATA" > "$METADATA"
fi

# ---------------------------
# 5) Downstream 컬럼
# ---------------------------
IFS=',' read -ra GROUP_COLUMNS <<< "${DOWN_STREAM_COLUMNS:-group}"

# ---------------------------
# 6) PIPELINE 실행
# ---------------------------
PIPELINE "$MANIFEST" "$OUTPUT_NAME" "$SEQ_MODE" "${GROUP_COLUMNS[@]}"

