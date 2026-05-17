#!/bin/bash
export PATH="$HOME/Master_thesis/bin/meme-5.5.9/bin:$PATH"

INPUT_DIR="../output/5_2_Files_for_motif_enrichment"
MEME_LIB="../data/5_HOCOMOCO_motifs_db_v5.meme"
GLOBAL_OUT="../output/5_3A_Global_motif_enrichment"

mkdir -p "$GLOBAL_OUT"

ALL_PS="$GLOBAL_OUT/ALL_PS_sequences.fasta"
ALL_NPS="$GLOBAL_OUT/ALL_NPS_sequences.fasta"
ALL_PS_DEDUP="$GLOBAL_OUT/ALL_PS_sequences_dedup.fasta"
ALL_NPS_DEDUP="$GLOBAL_OUT/ALL_NPS_sequences_dedup.fasta"

# Pool all PS and NPS sequences
cat "$INPUT_DIR"/*/PS_sequences.fasta > "$ALL_PS"
cat "$INPUT_DIR"/*/NPS_sequences.fasta > "$ALL_NPS"

# Deduplicate by header
awk '/^>/{if (seen[$0]++) skip=1; else skip=0} !skip' "$ALL_PS" > "$ALL_PS_DEDUP"
awk '/^>/{if (seen[$0]++) skip=1; else skip=0} !skip' "$ALL_NPS" > "$ALL_NPS_DEDUP"

# Report counts
PS_RAW=$(grep -c ">" "$ALL_PS")
NPS_RAW=$(grep -c ">" "$ALL_NPS")
PS_DEDUP=$(grep -c ">" "$ALL_PS_DEDUP")
NPS_DEDUP=$(grep -c ">" "$ALL_NPS_DEDUP")
echo "PS:  $PS_RAW → $PS_DEDUP (removed $((PS_RAW - PS_DEDUP)) duplicates)"
echo "NPS: $NPS_RAW → $NPS_DEDUP (removed $((NPS_RAW - NPS_DEDUP)) duplicates)"

# Run AME globally on deduplicated files
echo "Running global AME enrichment..."
ame --oc "$GLOBAL_OUT" \
    --control "$ALL_NPS_DEDUP" \
    "$ALL_PS_DEDUP" \
    "$MEME_LIB"

echo "Done. Results in $GLOBAL_OUT"