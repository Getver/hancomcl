#!/bin/bash
# microbiome shotgun sequencing read-based annotation PIPELINE
# Usage:
# bash RUN_PIPELINE.sh --fastq-dir <dir> --work-dir <dir> --name <name>
# bash RUN_PIPELINE.sh --fastq-dir <dir> --work-dir <dir> --name <name> --meta-data <file> --down-stream Columns1,Columns2 --service QC

set -euo pipefail

# ---------------------------
# 0) 기본값 초기화
# ---------------------------
FASTQ_DIR=""
WORK_DIR=""
KEY=""
META_DATA=""
SERVICE=""
DOWN_STREAM=false
DOWN_STREAM_COLUMNS=""

# ---------------------------
# 1) 옵션 파싱
# ---------------------------
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --name) KEY="$2"; shift 2 ;;
        --meta-data) META_DATA="$2"; shift 2 ;;
        --down-stream)
            DOWN_STREAM=true
            # 다음 인자가 존재하고 --로 시작하지 않으면 컬럼명으로 인식
            if [[ -n "${2:-}" && ! "$2" =~ ^-- ]]; then
                DOWN_STREAM_COLUMNS="$2"
                shift 2   # 옵션과 컬럼명 둘 다 shift
            else
                shift
            fi
            ;;
        --service) SERVICE="$2"; shift 2 ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: bash RUN_PIPELINE.sh --fastq-dir <dir> --work-dir <dir> --name <name> [--meta-data <file>] [--down-stream Columns1,Columns2] [--service <name>]"
            exit 1 ;;
    esac
done

# ---------------------------
# 2) 필수값 체크
# ---------------------------
if [[ -z "$FASTQ_DIR" || -z "$WORK_DIR" || -z "$KEY" ]]; then
    echo "[ERROR] --fastq-dir, --work-dir, and --name must all be specified."
    exit 1
fi

# ---------------------------
# 3) 출력 폴더 생성
# ---------------------------
OUTPUT_PATH="$WORK_DIR/$KEY"
mkdir -p "$OUTPUT_PATH"

# mv $META_DATA $KEY/

LOG_DIR="$OUTPUT_PATH/logs"
mkdir -p "$LOG_DIR"
ORDER_FILE="$LOG_DIR/order.tsv"


info_msg="FASTQ_DIR: $FASTQ_DIR
WORK_DIR: $WORK_DIR
OUTPUT_FOLDER: $KEY
META_DATA: ${META_DATA:-none}
SERVICE: ${SERVICE:-none}
DOWN_STREAM: $DOWN_STREAM
DOWN_STREAM_COLUMNS: ${DOWN_STREAM_COLUMNS:-none}"

echo "$info_msg" >> "$ORDER_FILE"


# ---------------------------
# 4) Docker 실행 명령 조립
# ---------------------------
DOCKER_CMD=(
    docker run --rm
    -v "$OUTPUT_PATH:/home/data"
    -v "$FASTQ_DIR:/work/raw_data"
    -v "/disk0/sm/microbiome/16S/docker/DB:/work/DB"
    -v "/disk0/sm/reference/ncbi:/work/ncbi"
)

# ---------------------------
# 5) META_DATA 처리
# ---------------------------
TARGET_META="$OUTPUT_PATH/file_list.txt"

if [[ -n "$META_DATA" ]]; then
    if [[ "$META_DATA" != "$TARGET_META" ]]; then
        cp "$META_DATA" "$TARGET_META"
        echo "[INFO] META_DATA 복사됨 → $TARGET_META"
    else
        echo "[INFO] META_DATA 이미 작업폴더에 존재 → $TARGET_META"
    fi
    DOCKER_CMD+=( -e "RAW_METADATA=/home/data/$(basename "$TARGET_META")" )
else
    echo "[INFO] META_DATA 없음 (옵션 미제공)"
fi

# ---------------------------
# 6) 기타 옵션 전달
# ---------------------------
[[ "$DOWN_STREAM" == true ]] && DOCKER_CMD+=( -e "DOWN_STREAM=true" )
[[ -n "$DOWN_STREAM_COLUMNS" ]] && DOCKER_CMD+=( -e "DOWN_STREAM_COLUMNS=$DOWN_STREAM_COLUMNS" )
[[ -n "$SERVICE" ]] && DOCKER_CMD+=( -e "SERVICE=$SERVICE" )

# ---------------------------
# 7) 이미지와 인자 추가
# ---------------------------
DOCKER_CMD+=( "sm:qiime_0.2" "$KEY" )

# ---------------------------
# 8) 실행
# ---------------------------
echo "[INFO] 실행 명령:"
echo "${DOCKER_CMD[*]}"
"${DOCKER_CMD[@]}"


