#!/usr/bin/env python3
"""
Batch GC% extraction pipeline.

What this script does (per TF, per cross-validation fold):
  1) Find test BED files named like 'test_FXX.bed' under '<TF>/crossval/fold_XX/'.
  2) Run 'bedtools nuc' against hg19 to compute nucleotide composition.
  3) Parse the 'bedtools nuc' output and write a compact TSV with: chrom, start, end, GC%.

Inputs/assumptions:
  - Project layout: ~/msc_project/data/<TF>/crossval/fold_XX/test_FXX.bed
  - Reference FASTA at: ~/genomes/hg19.fa   (NOTE: 'hg19.fa.fai' should exist alongside.)
  - Tools on PATH: 'bedtools'
  - 'bedtools nuc' output format is expected such that the GC fraction is in column 5

Outputs (per fold):
  - '<fold_dir>/bg_BXX.nuc.txt' : raw 'bedtools nuc' output
  - '<fold_dir>/bg_BXX.gc.txt'  : tab-seperated: chrom, start, end, GC%

Notes:
  - If the 'bedtools nuc' output format differs across versions, adjust the index used for GC
  - Folders without the expected structure/files are skipped
"""
from pathlib import Path
import subprocess
import sys

DATA_ROOT = Path.home() / "msc_project" / "data" # Path to the data directory
CRG_ROOT = Path.home() / "genomes" # Path to the CRG genome directory

# Loop over each TF folder under data_root
for tfdir in sorted(DATA_ROOT.iterdir()):
    crossval = tfdir / "crossval"
    # skip anything that doesn't have a cross-validation subdirectory
    if not crossval.is_dir():
        continue

    print(f"==== Processing TF: {tfdir.name} ====")
    # iterate over each fold directory
    for fold_dir in sorted(crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue

        # find test_FXX.bed in this fold
        for bed_path in fold_dir.glob("test_F[0-9][0-9].bed"):
            print(f"Processing {bed_path.name} in {fold_dir.name}")
            
            # derive a background name suffix from the fold number
            fold_number = fold_dir.name.split('_')[-1]
            name = f"bg_B{fold_number.zfill(2)}"  # e.g. "fold_01" → "bg_B01"
            
            # run bedtools nuc to compute nucleotide composition for each interval
            print("Running bedtools nuc for GC% calculation...")
            cmd = f"bedtools nuc -fi {CRG_ROOT / 'hg19.fa'} -bed {bed_path} > {fold_dir}/{name}.nuc.txt"
            subprocess.run(cmd, shell=True, check=True)
            print(f"==== Bedtools nuc completed for {bed_path.name} in {fold_dir.name} ====")
            
            # Extract chr, start, end, GC% from nuc output
            print(f"Extracting chr, start, end, GC% from {name}.nuc.txt...")
            print(f"Writing to {fold_dir / f'{name}.gc.txt'}")
            
            # open the raw nuc output for reading and the GC summary file for writing 
            with open(fold_dir / f"{name}.nuc.txt") as f, open(fold_dir / f"{name}.gc.txt", "w") as out:
                for line in f:
                    # skip commented header lines
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


            