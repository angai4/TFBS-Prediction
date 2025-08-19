#!/usr/bin/env python3
"""
Compute GC% for FNN training windows.

For each TF directory under ~/msc_project/data that contains `fnn_crossval/`,
this script:

  1) Finds each fold folder:         fnn_crossval/fold_XX/
  2) Looks for a BED file named:     train_pos_TXX.bed
  3) Runs `bedtools nuc` against hg19 to compute per-window nucleotide stats.
  4) Writes the raw output to:       train_pos_TXX.nuc.txt
  5) Extracts chr, start, end, GC%   -> train_pos_TXX.gc.txt   (GC as integer %)

Requirements:
  - bedtools installed and on PATH
  - Reference genome at ~/genomes/hg19.fa (ideally indexed with `samtools faidx`)

Notes:
  - `bedtools nuc` emits several columns; we use:
        0: chrom, 1: start, 2: end, 4: GC fraction (0..1)
    Lines beginning with '#' are headers and are skipped.
"""

from pathlib import Path
import subprocess
import sys

DATA_ROOT = Path.home() / "msc_project" / "data" # Path to the data directory
CRG_ROOT = Path.home() / "genomes" # Path to the CRG genome directory

# Loop over each TF folder under data_root
for tfdir in sorted(DATA_ROOT.iterdir()):
    fnn_crossval = tfdir / "fnn_crossval"
    if not fnn_crossval.is_dir():
        continue

    print(f"==== Processing TF: {tfdir.name} ====")
    for fold_dir in sorted(fnn_crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue

        # find train_pos_TXX.bed in this fold
        for bed_path in fold_dir.glob("train_pos_T[0-9][0-9].bed"):
            print(f"Processing {bed_path.name} in {fold_dir.name}")
            fold_number = fold_dir.name.split('_')[-1]
            name = f"T{fold_number.zfill(2)}"  # e.g. "fold_01" → "T01"
            print("Running bedtools nuc for GC% calculation...")
            cmd = f"bedtools nuc -fi {CRG_ROOT / 'hg19.fa'} -bed {bed_path} > {fold_dir}/train_pos_{name}.nuc.txt"
            subprocess.run(cmd, shell=True, check=True)
            print(f"==== Bedtools nuc completed for {bed_path.name} in {fold_dir.name} ====")
            # Extract chr, start, end, GC% from nuc output
            print(f"Extracting chr, start, end, GC% from train_pos_{name}.nuc.txt...")
            print(f"Writing to {fold_dir / f'train_pos_{name}.gc.txt'}")
            with open(fold_dir / f"train_pos_{name}.nuc.txt") as f, open(fold_dir / f"train_pos_{name}.gc.txt", "w") as out:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    if len(parts) < 12:
                        print(f"  WARNING: unexpected format in {line.strip()}", file=sys.stderr)
                        continue
                    chrom, start, end, gc = parts[0], parts[1], parts[2], int(round(float(parts[4]) * 100))
                    out.write(f"{chrom}\t{start}\t{end}\t{gc}\n")
            print(f"==== GC% extraction completed for {name} ====")

print("All TFs processed.")


            