Ncbi_blast () {

    export BLASTDB=/work/ncbi/16S_ribosomal_RNA/

    blastn \
    -query filtered/seq/denoise_noMT_seq.fasta \
    -db 16S_ribosomal_RNA \
    -out taxonomy/blast_results.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads 50

}
