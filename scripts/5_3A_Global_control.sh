OUTPUT_DIR="../output/5_3A_Global_motif_enrichment"
MEME_LIB="../data/5_HOCOMOCO_motifs_db_v5.meme"
ALL_PS="$OUTPUT_DIR/ALL_PS_sequences_dedup.fasta"
SHUFFLE_OUT="$OUTPUT_DIR/control"

mkdir -p "$SHUFFLE_OUT"

echo "Running AME with shuffle control..."
ame --oc "$SHUFFLE_OUT" \
    --control --shuffle-- \
    "$ALL_PS" \
    "$MEME_LIB"

echo "Done. Results in $SHUFFLE_OUT"