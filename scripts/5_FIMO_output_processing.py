import os
import pandas as pd
from glob import glob
from Bio import SeqIO
from pathlib import Path

fasta_path = '../output/1_2_FASTA_per_enh'
motifs_mapped = '../output/2_1_Motifs_per_enhancer'
files = sorted(glob(f"{motifs_mapped}/*_motifs.csv"))


for file in files:
    name = Path(file).name.replace("_motifs.csv", "")
    df = pd.read_csv(file, sep=';', header=0)
    
    #motif length
    df['motif_length'] = abs(df['stop'] - df['start']) + 1
    
    #enhancer length
    fasta = f'{fasta_path}/{name}.fa'
    first = next(SeqIO.parse(fasta, "fasta"), None)
    df["enh_length"] = len(first.seq)
    df['enh_name'] = name
    
    df.to_csv(f'{motifs_mapped}/{name}.csv', index=False, sep=';')