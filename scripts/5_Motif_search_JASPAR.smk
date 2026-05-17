#TO RUN
#snakemake -s 2_1_Motif_search.smk -j 16 -n -p #dry-run
#snakemake -s 2_1_Motif_search.smk -j 16 &> ../logs/2_1_Motif_search.log

#libraries
import os
import pandas as pd
import glob
from pymemesuite.fimo import FIMO 
from pymemesuite.common import MotifFile 
from pathlib import Path
from pymsaviz import MsaViz
from Bio import AlignIO, Align
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import re
import glob

#data
meme_lib_url = "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
ref = 'hg38'

# list of enhancers with available trees (and MSAs), formed by rule list_tree_enhancers in the 3d snakemake file
fasta_dir = "../output/1_2_FASTA_per_enh"
aligned_enh = sorted(
    os.path.splitext(os.path.basename(p))[0]
    for p in glob.glob(os.path.join(fasta_dir, "*.fa"))
)

JASPAR_URL = "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"

#rules
rule all:
    input:
        expand("../output/2_1_Motifs_per_enhancer/{enhancer}_motifs.csv", enhancer=aligned_enh)

rule download_meme_library:
    """Download JASPAR (MEME format) as a reproducible step."""
    output:
        "../data/5_JASPAR_motifs_database.txt"
    shell:
        r"""
        mkdir -p $(dirname {output})
        curl -L "{JASPAR_URL}" -o {output}
        """

rule motif_search:
    """ 
    With available database of vertebrates motifs, the rule look through fasta sequences before multiple sequence alignment
    and identify all motifs. As an output, it provides .csv file per enhancer which lists all motifs in the sequences
    in which seqience it was found, statistical scoring, etc.
    """  
    input:
        motifs="../data/2_JASPAR_motifs_database.txt",
        fasta="../output/1_2_FASTA_per_enh/{enhancer}.fa"
    output:
        "../output/2_1_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt"
    params:
        ref_latin = ref,
        base_dir = "../output/2_1_Motifs_per_enhancer/uncurated",
        thresh = "1e-4",
        markov_order = 1
    shell:
        r"""
        mkdir -p {params.base_dir}

        outdir="{params.base_dir}/{wildcards.enhancer}"
        mkdir -p "$outdir"

        ref_fa="$outdir/{wildcards.enhancer}.ref.fa"
        bg="$outdir/{wildcards.enhancer}.bg"

        # 1) Extract only the reference organism sequence
        awk -v ref="{params.ref_latin}" '
            BEGIN {{p=0}}
            /^>/ {{
                p = (index($0, ref) > 0)
            }}
            p {{print}}
        ' {input.fasta} > "$ref_fa"

        # 2) Build background model (Markov)
        # fasta-get-markov writes background file; keep it in per-enhancer outdir
        fasta-get-markov -m {params.markov_order} "$ref_fa" "$bg" >/dev/null 2>&1

        # 3) adjusting p-value threshold
        ref_len=$(
        awk '
            /^>/ {{next}}
            {{ gsub(/[ \t\r\n]/, "", $0); len += length($0) }}
            END {{ print len+0 }}
        ' "$ref_fa"
        )

        p0="{params.thresh}"   # e.g. 1e-4
        p_eff=$(
        awk -v L="$ref_len" -v p="$p0" '
            BEGIN {{
                if (L > 1000) printf "%.10g", (p * 1000.0 / L);
                else printf "%.10g", p;
            }}
        '
        ) #recalculation of the p-value threshold

        # 4) Run FIMO; do not fail the rule if it fails
        fimo --oc "$outdir" --thresh $p_eff --bgfile "$bg" {input.motifs} "$ref_fa" >/dev/null 2>&1

        # Pick the produced FIMO table (format may vary by version)
        src=""
        if [ -f "$outdir/fimo.tsv" ]; then
          src="$outdir/fimo.tsv"
        elif [ -f "$outdir/fimo.txt" ]; then
          src="$outdir/fimo.txt"
        fi

        # If no output produced -> create empty file, but do not crash
        if [ -z "$src" ] || [ ! -s "$src" ]; then
          echo "FIMO produced no result table for {wildcards.enhancer}"
          : > {output}
          exit 0
        fi

        mv "$src" {output}
        """

rule curating_meme_output:
    """ 
    Initial output of fimo doesn't provide normal naming of found motifs. The dataset contains only motifs IDs which
    is not really convenient in analysis proceeding. Thus, the rule takes JASPAR database, where each motif is named as
    "MOTIF MA0069.1 PAX6", dissect IDs and proper name, goes through all output csvs from the previous rule and replace
    IDs by proper name. 
    """  
    input:
        motifs_txt = "../output/2_1_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt",
        meme_lib   = "../data/2_JASPAR_motifs_database.txt"
    output:
        motifs_csv = "../output/2_1_Motifs_per_enhancer/{enhancer}_motifs.csv"
    run:
        # dictionary "MAxxxx.x -> motif_name"
        mapping = {}
        with open(input.meme_lib, encoding="utf-8") as fh:
            for line in fh:
                m = re.match(r'^MOTIF\s+(\S+)(?:\s+(.+))?$', line.strip())
                if m:
                    motif_id = m.group(1)
                    alt = (m.group(2) or "").strip()
                    if alt:
                        mapping[motif_id] = alt

        # reading FIMO output
        df = pd.read_csv(input.motifs_txt, sep=None, engine="python")

        rename_map = {
            "#pattern name": "motif_id",
            "motif_id": "motif_id",           
            "sequence name": "organism",
            "matched sequence": "sequence",
            "p-value": "p_value",
            "q-value": "q_value",
            "start": "start",
            "stop": "stop",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        # column with readble motif name
        if "motif_id" in df.columns:
            df["motif_name"] = df["motif_id"].map(mapping).fillna(df["motif_id"])
        else:
            df["motif_name"] = ""

        df.to_csv(output.motifs_csv, sep=';', index=False)

