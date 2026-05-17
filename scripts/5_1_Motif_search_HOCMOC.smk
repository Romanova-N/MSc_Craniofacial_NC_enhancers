# snakemake -s 2_1H_Motif_search_HOCMOC.smk -j 16 -n -p  # dry-run
# snakemake -s 2_1H_Motif_search_HOCMOC.smk -j 16 &> ../logs/2_1H_Motif_search.log

import os
import re
import glob
import pandas as pd

ref = 'hg38'

fasta_dir = "../output/1_2_FASTA_per_enh"
aligned_enh = sorted(
    os.path.splitext(os.path.basename(p))[0]
    for p in glob.glob(os.path.join(fasta_dir, "*.fa"))
)

rule all:
    input:
        expand("../output/2_1H_Motifs_per_enhancer/{enhancer}_motifs.csv", enhancer=aligned_enh)


rule motif_search:
    """
    Run FIMO against the HOCOMOCO MEME database on the hg38 reference
    sequence for each enhancer. Threshold is length-adjusted. Produces
    a raw fimo.tsv (or an empty file if FIMO finds nothing).
    """
    input:
        motifs = "../data/5_HOCOMOCO_motifs_db.meme",
        fasta  = "../output/1_2_FASTA_per_enh/{enhancer}.fa"
    output:
        "../output/2_1H_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt"
    params:
        ref_tag      = ref,
        base_dir     = "../output/2_1H_Motifs_per_enhancer/uncurated",
        thresh       = "1e-4",
        markov_order = 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.base_dir}

        outdir="{params.base_dir}/{wildcards.enhancer}"
        mkdir -p "$outdir"

        ref_fa="$outdir/{wildcards.enhancer}.ref.fa"
        bg="$outdir/{wildcards.enhancer}.bg"
        fimo_log="$outdir/{wildcards.enhancer}.fimo.log"

        # 1) Extract only the reference organism sequence
        awk -v ref="{params.ref_tag}" '
            BEGIN {{ p=0 }}
            /^>/ {{ p = (index($0, ref) > 0) }}
            p {{ print }}
        ' {input.fasta} > "$ref_fa"

        # 2) Build background model
        fasta-get-markov -m {params.markov_order} "$ref_fa" "$bg" \
            > "$outdir/fasta-get-markov.log" 2>&1

        # 3) Length-adjusted p-value threshold
        ref_len=$(awk '
            /^>/ {{ next }}
            {{ gsub(/[ \t\r\n]/, "", $0); len += length($0) }}
            END {{ print len+0 }}
        ' "$ref_fa")

        p_eff=$(awk -v L="$ref_len" -v p="{params.thresh}" '
            BEGIN {{
                if (L > 1000) printf "%.10g\n", (p * 1000.0 / L);
                else printf "%.10g\n", p;
            }}
        ')

        # 4) Run FIMO; stderr goes to per-enhancer log
        fimo --oc "$outdir" --thresh "$p_eff" --bgfile "$bg" \
            {input.motifs} "$ref_fa" > "$fimo_log" 2>&1 || true

        # 5) Collect output
        src=""
        [ -f "$outdir/fimo.tsv" ] && [ -s "$outdir/fimo.tsv" ] && src="$outdir/fimo.tsv"
        [ -z "$src" ] && [ -f "$outdir/fimo.txt" ] && [ -s "$outdir/fimo.txt" ] && src="$outdir/fimo.txt"

        if [ -z "$src" ]; then
            echo "FIMO produced no hits for {wildcards.enhancer}" >&2
            : > {output}
            exit 0
        fi

        mv "$src" {output}
        """


rule curating_meme_output:
    """
    HOCOMOCO (human-only DB) MEME headers have the form:
        MOTIF AHR.H14CORE.0.P.B
    The TF symbol is everything before the first '.' — no _HUMAN suffix,
    no second alt-name field. FIMO puts the full model ID in motif_id.
    We extract the symbol and keep the full ID for traceability.
    Guards against empty FIMO output files.
    """
    input:
        motifs_txt = "../output/2_1H_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt"
    output:
        motifs_csv = "../output/2_1H_Motifs_per_enhancer/{enhancer}_motifs.csv"
    run:
        # Handle empty FIMO output gracefully
        output_columns = ['motif_id', 'organism', 'start', 'stop', 'strand', 'score', 'p_value', 'q_value', 'sequence']

        if os.path.getsize(input.motifs_txt) == 0:
            pd.DataFrame(columns=output_columns).to_csv(output.motifs_csv, sep=';', index=False)
            return

        # Read FIMO TSV; skip comment lines starting with '#'
        df = pd.read_csv(
            input.motifs_txt,
            sep='\t',
            engine='python'
        )
        
        if df.empty:
            pd.DataFrame(columns=output_columns).to_csv(output.motifs_csv, sep=';', index=False)
            return

        df = df.rename(columns={
                '#pattern name':  'motif_id',
                'sequence name':  'organism',
                'matched sequence': 'sequence',
                'p-value':        'p_value',
                'q-value':        'q_value',
            })
        
        # Extract TF symbol from HOCOMOCO model ID:
        # 'AHR.H14CORE.0.P.B' in motif_name = 'AHR'
        if "motif_id" not in df.columns:
            raise ValueError(f'There is no motif_id column in {input.motifs_txt}. '
                             f'Found: {list(df.columns)}')
        df["motif_name"] = df["motif_id"].str.split('.').str[0]

        df.to_csv(output.motifs_csv, sep=';', index=False)
