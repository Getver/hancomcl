# source /disk0/sm/microbiome/docker/script/*
source /work/module/ONESAMPLE_PIPELINE.sh
source /work/module/QC_fastp.sh
source /work/module/HostRemove_bowtie2.sh
source /work/module/Preprocessing_FASTQ.sh
source /work/module/Taxonomy_metaphlan.sh
source /work/module/Functional_humann.sh
source /work/module/Taxonomy_kraken.sh
source /work/module/Abundance_bracken.sh
source /work/module/Result_merge.sh
source /work/module/Diversity_AB.sh
source /work/module/Tree_krona.sh

set -euo pipefail

## database 폴더 업데이트(인식시켜주기) : 최초 실행시에만 필요
humann_config --update database_folders utility_mapping /work/DB/humann_db/utility_mapping && apt-get install -y libglpk40 glpk-utils

if [[ -n "$READ1" && -n "$READ2" ]]; then
    # 절대경로가 들어와도 파일명만 추출
    R1=$(basename "$READ1")
    R2=$(basename "$READ2")

    # SAMPLE 추출
    if [[ "$R1" =~ _R1_ ]]; then
        SAMPLE="${R1%%_R1_*}"
    elif [[ "$R1" =~ _R1\. ]]; then
        SAMPLE="${R1%%_R1.*}"
    elif [[ "$R1" =~ _1\.fastq$ ]]; then
        SAMPLE="${R1%%_1.fastq*}"
    elif [[ "$R1" =~ _1\.fastq\.gz$ ]]; then
        SAMPLE="${R1%%_1.fastq.gz}"
    else
        SAMPLE="unknown"
    fi

    echo "[INFO] READ1 / READ2 입력됨 - $R1, $R2"
    echo "[INFO] SAMPLE: $SAMPLE"

    ONESAMPLE_PIPELINE "$SAMPLE" "/work/raw_data/${R1}" "/work/raw_data/${R2}"
else
    echo "[INFO] READ1/READ2가 없어 자동 매칭 수행 중..."

    FILE_LIST=$(ls /work/raw_data/*fastq* | xargs -n 1 basename)

    for R1 in $(echo "$FILE_LIST" | grep -E '_R?1(_|[._])|_1\.fastq'); do
        if [[ "$R1" =~ _R1_ ]]; then
            R2="${R1/_R1_/_R2_}"
            SAMPLE="${R1%%_R1_*}"
        elif [[ "$R1" =~ _R1\. ]]; then
            R2="${R1/_R1./_R2.}"
            SAMPLE="${R1%%_R1.*}"
        elif [[ "$R1" =~ _1\.fastq$ ]]; then
            R2="${R1/_1.fastq/_2.fastq}"
            SAMPLE="${R1%%_1.fastq}"
        elif [[ "$R1" =~ _1\.fastq\.gz$ ]]; then
            R2="${R1/_1.fastq.gz/_2.fastq.gz}"
            SAMPLE="${R1%%_1.fastq.gz}"
        else
            echo "[WARN] Unknown R1 format: $R1"
            continue
        fi

        if echo "$FILE_LIST" | grep -q "$R2"; then
            echo "샘플 자동 매칭:"
            echo "   SAMPLE: $SAMPLE"
            echo "   R1: $R1"
            echo "   R2: $R2"

            ONESAMPLE_PIPELINE ${SAMPLE} /work/raw_data/${R1} /work/raw_data/${R2}

        else
            echo "[WARN] 자동 매칭 실패 샘플: $R1 (예상: $R2)"
        fi
    done
    Result_merge
    Diversity_AB
fi
