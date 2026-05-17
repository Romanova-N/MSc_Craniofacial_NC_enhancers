import sys
import glob
import os
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
from bx.align import maf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# output argument
if len(sys.argv) > 4:
    output_dir = sys.argv[4]
    os.makedirs(output_dir, exist_ok=True)
else:
    output_dir = sys.argv[3]

#functions
def cure_header(src: str) -> str:
    return src.split(".", 1)[0]


def maf_2_fasta(
    maf_path,
    species: set,
    concat=True,
    keep_gaps=True,
    fill_missing="-",
    out_fasta=None
):
    """
    Convert a (subset) MAF alignment into a FASTA file with one sequence per species.

    This version is robust to a common bigMaf/MAF quirk:
    - The SAME species may appear multiple times within the SAME MAF block.
      Naively appending all occurrences makes that species artificially longer.
      Here we choose ONE "best" sequence per species per block, then append exactly once.

    "Best" duplicate selection strategy (per block, per species):
    - Prefer the component with the FEWEST gap characters ('-')
      (equivalently: the MOST non-gap aligned bases).

    Also:
    - Species missing in a block are padded with `fill_missing * block_len`
      so all output sequences keep the same alignment length after concatenation.

    Returns:
      (out_fasta_path, n_records_written)
    """
    maf_path = Path(maf_path)

    if out_fasta is None:
        out_fasta = os.path.join(
            output_dir,
            os.path.basename(str(maf_path)).replace(".subset.maf", ".fa")
        )

    # Accumulate per-block strings for each species
    seq_per_sp = defaultdict(list)

    with maf_path.open("r") as fh:
        reader = maf.Reader(fh)

        for block in reader:
            block_len = None

            # For this block, keep only ONE sequence per species (resolve duplicates)
            block_best = {}  # sp -> best aligned string

            for comp in block.components:
                if comp.text is None:
                    continue  # skip empty components

                if block_len is None:
                    block_len = len(comp.text)

                sp = cure_header(comp.src)
                seq = comp.text

                # If species not requested, skip (optional but can save memory/time)
                if sp not in species:
                    continue

                # First occurrence
                if sp not in block_best:
                    block_best[sp] = seq
                    continue

                # Duplicate occurrence: pick the "best" one (fewest gaps)
                # (tie-breaker: keep the current best)
                if seq.count("-") < block_best[sp].count("-"):
                    block_best[sp] = seq

            if block_len is None:
                continue  # block had no sequences

            # Append exactly once per present species
            present_sp = set(block_best.keys())
            for sp, seq in block_best.items():
                seq_per_sp[sp].append(seq)

            # Pad missing species so alignment lengths match after concatenation
            missing_sp = species - present_sp
            for sp in missing_sp:
                seq_per_sp[sp].append(fill_missing * block_len)

    # Build SeqRecords for all requested species
    records = []
    for sp in species:
        seq_str = "".join(seq_per_sp[sp])

        # Optionally remove gaps to output ungapped sequences
        if not keep_gaps:
            seq_str = seq_str.replace("-", "")

        records.append(SeqRecord(Seq(seq_str), id=sp, description=""))

    SeqIO.write(records, out_fasta, "fasta")
    return out_fasta, len(records)

#workflow
ref_sp = sys.argv[1]

with open(sys.argv[2], "r") as f:
    species_list = [line.strip() for line in f if line.strip()] # list of species ids
    species_list.append(ref_sp)

input_path = sys.argv[3]

if os.path.isdir(input_path):
    maf_directory = input_path
    maf_files = glob.glob(f"{maf_directory}/*.subset.maf")

    for maf_file in tqdm(maf_files):
        name, nseq = maf_2_fasta(maf_file, set(species_list))
        print(f"[DONE] {name}, {nseq} sequences in the file")

elif os.path.isfile(input_path):
    name, nseq = maf_2_fasta(input_path, set(species_list))
    print(f"[DONE] {name}, {nseq} sequences in the file")

else:
    print("[ERROR] Check your arguments!")
