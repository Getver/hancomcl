sanitize_metadata() {
    awk '
    BEGIN { FS="\t"; OFS="\t" }
    NR==1 { print; next }

    {
        for (i=1; i<=NF; i++) {
            val = $i

            gsub(/[\/]/, "-", val)   # / → -
            gsub(/[ ]+/, "_", val)   # space → _
            gsub(/[^A-Za-z0-9_\-\.]/, "", val)  # remove other weird chars

            $i = val
        }
        print
    }
    ' "$1"
}