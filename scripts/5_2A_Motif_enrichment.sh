INPUT_DIR="../output/5_2_Files_for_motif_enrichment"
OUTPUT_DIR="../output/5_2A_Motif_enrichment"
MEME_LIB="../data/5_HOCOMOCO_motifs_db_v5.meme"
GENE_LIST="../data/5_2_target_genes.txt"
GENES=$(cat "$GENE_LIST")

mkdir -p "$OUTPUT_DIR"

for gene in $GENES; do
    GENE_OUT="$OUTPUT_DIR/$gene"
    mkdir -p "$GENE_OUT"

    PS_FASTA="$INPUT_DIR/$gene/PS_sequences.fasta"
    NPS_FASTA="$INPUT_DIR/$gene/NPS_sequences.fasta"

    if [[ ! -s "$PS_FASTA" || ! -s "$NPS_FASTA" ]]; then
        echo "WARNING: $gene — missing or empty FASTA, skipping." >&2
        continue
    fi

    echo "Running AME enrichment for $gene..."
    ame --oc "$GENE_OUT" --control "$NPS_FASTA" "$PS_FASTA" "$MEME_LIB"
done

echo "Finished"