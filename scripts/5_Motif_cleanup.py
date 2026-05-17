import pandas as pd
import glob
import os
import numpy as np
from pathlib import Path
from Bio import SeqIO
from itertools import combinations

motifs_mapped = '../output/2_1H_Motifs_per_enhancer'
fasta_path = '../output/1_2_FASTA_per_enh'
files = glob.glob(f"{motifs_mapped}/*.csv")

# list of TFs expressed in our dataset
with open('../data/5_targeted_TFs.txt', 'r') as handle:
    TFs_lib = [line.strip() for line in handle if line.strip()]

output_dir = '../output/2_2H_Motifs_per_enhancer'
os.makedirs(output_dir, exist_ok=True)

# functions
def filter_tf_hits_by_overlap_single_enhancer(
    df,tf_col="motif_name",
    start_col="start",stop_col="stop",
    len_col="motif_length",score_col="score",
    qval_col="q_value",min_recip_overlap=0.7):

    '''
    The function filter situations when two motifs fall into the same or almost same region
    Only one motif will be saved based on either FIMO score (if the same motif) or
    q-value (if the motifs are different)
    '''

    work = df.copy().reset_index(drop=False).rename(columns={"index": "_orig_idx"})
    to_remove = set()

    for i, j in combinations(work.index, 2):
        r1 = work.loc[i]
        r2 = work.loc[j]

        if r1["_orig_idx"] in to_remove or r2["_orig_idx"] in to_remove:
            continue

        overlap = max(0, min(r1[stop_col], r2[stop_col]) - max(r1[start_col], r2[start_col]) + 1)
        recip1 = overlap / r1[len_col]
        recip2 = overlap / r2[len_col]

        if recip1 < min_recip_overlap or recip2 < min_recip_overlap:
            continue

        same_tf = r1[tf_col] == r2[tf_col]

        if same_tf:
            # same motif -> higher score wins
            if r1[score_col] > r2[score_col]:
                drop_idx = r2["_orig_idx"]
            elif r2[score_col] > r1[score_col]:
                drop_idx = r1["_orig_idx"]
            else:
                # equal scores: keep both
                continue

        else:
            # different motifs -> lower q-value wins
            if r1[qval_col] < r2[qval_col]:
                drop_idx = r2["_orig_idx"]
            elif r2[qval_col] < r1[qval_col]:
                drop_idx = r1["_orig_idx"]
            else:
                # equal q-values: keep both
                continue

        to_remove.add(drop_idx)

    return work.loc[~work["_orig_idx"].isin(to_remove)].drop(columns=["_orig_idx"])

q_threshold = 0.05

for file in files:
    name = Path(file).name.replace("_motifs.csv", "")
    df = pd.read_csv(file, sep=';', header=0)
    
    #motif length
    df['motif_length'] = abs(df['stop'].astype(int) - df['start'].astype(int)) + 1
    
    #enhancer length
    fasta = f'{fasta_path}/{name}.fa'
    first = next(SeqIO.parse(fasta, "fasta"), None)
    df["enh_length"] = len(first.seq)
    df['enh_name'] = name
    
    # keep only uppercase motif names (so human-derived)
    df = df[df["motif_name"] == df["motif_name"].str.upper()].copy()

    # keep only TFs expressed in our dataset
    df = df[df['motif_name'].isin(TFs_lib)].copy()

    # filter by q-value
    df = df[df['q_value'] < q_threshold].copy()

    # filter hits falling in the same region
    df_to_save = filter_tf_hits_by_overlap_single_enhancer(df)

    df_to_save.to_csv(f'{output_dir}/{name}.csv', sep=';', index=False)