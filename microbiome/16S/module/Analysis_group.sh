Analysis_group() {
    local GROUP_COLUMNS=("$@")

    # Alpha diversity
    mkdir -p diversity/alpha_group
    alpha_metrics=("faith_pd_vector" "evenness_vector" "shannon_vector")
    for metric in "${alpha_metrics[@]}"; do
        IN_QZA="diversity/${metric}.qza"
        OUT_QZV="diversity/alpha_group/${metric}-group-significance.qzv"
        HTML_DIR="diversity/alpha_group/${metric}-html"

        qiime diversity alpha-group-significance \
            --i-alpha-diversity "$IN_QZA" \
            --m-metadata-file "$METADATA" \
            --o-visualization "$OUT_QZV"

        mkdir -p "$HTML_DIR"
        tmpdir=$(mktemp -d)
        unzip -q -o "$OUT_QZV" -d "$tmpdir"
        uuid_dir=$(find "$tmpdir" -mindepth 1 -maxdepth 1 -type d | head -n1)
        [[ -d "$uuid_dir/data" ]] && mv "$uuid_dir/data/"* "$HTML_DIR/"
        rm -rf "$tmpdir"
    done

    # Beta diversity
    mkdir -p "diversity/beta_group"

    # find 결과를 null 문자로 처리하여 공백 문제 방지
    while IFS= read -r -d '' DIST_MATRIX; do
        base=$(basename "$DIST_MATRIX" _distance_matrix.qza)

        for col in "${GROUP_COLUMNS[@]}"; do
            safe_col=$(echo "$col" | tr ' ' '_' | tr -cd '[:alnum:]_-')
            OUT_QZV="diversity/beta_group/${base}-${safe_col}-significance.qzv"

            echo "[INFO] Processing $DIST_MATRIX → $OUT_QZV (column: $col)"
            qiime diversity beta-group-significance \
                --i-distance-matrix "$DIST_MATRIX" \
                --m-metadata-file "$METADATA" \
                --m-metadata-column "$col" \
                --p-method permanova \
                --o-visualization "$OUT_QZV"

            # HTML 추출
            HTML_DIR="diversity/beta_group/${base}-${safe_col}-significance-html"
            mkdir -p "$HTML_DIR"
            tmpdir=$(mktemp -d)
            unzip -q -o "$OUT_QZV" -d "$tmpdir"
            uuid_dir=$(find "$tmpdir" -mindepth 1 -maxdepth 1 -type d | head -n1)
            if [[ -d "$uuid_dir/data" ]]; then
                mv "$uuid_dir/data/"* "$HTML_DIR/" 2>/dev/null || true
            fi
            rm -rf "$tmpdir"
        done
    done < <(find "diversity" -type f -name "*_distance_matrix.qza" -print0)

    echo "[INFO] Analysis_group pipeline finished."
}
