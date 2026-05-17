URL=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.bigMaf
BED="../data/0_enhancer_coordinates.bed"
OUTDIR="../output/1_MAF_per_enh"
SPECIES="../data/0_target_species.txt"

mkdir -p "$OUTDIR"

sanitize() {
  echo "$1" | tr ' /' '__' | tr -cd 'A-Za-z0-9._-'
}

TOTAL=$(grep -cv '^[[:space:]]*$' "$BED")

pv -l -s "$TOTAL" "$BED" | \
while IFS=$'\t' read -r chr start end name || [[ -n "${chr:-}" ]]; do
  [[ -z "${chr:-}" ]] && continue
  name_s=$(sanitize "${name:-${chr}_${start}_${end}}")

  bigBedToBed "$URL" stdout -chrom="$chr" -start="$start" -end="$end" \
  | cut -f4 | tr ';' '\n' \
  | awk -f 1_1_maf_subset.awk "$SPECIES" - \
  > "${OUTDIR}/${name_s}.subset.maf"

done