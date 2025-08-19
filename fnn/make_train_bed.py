#!/usr/bin/env python3
"""
Create BED files for FNN training windows from FASTA headers.

For each TF under the provided data root, this scans:
    data/<TF>/fnn_crossval/fold_XX/train_pos_TXX.fa

and writes, alongside each FASTA, a BED file:
    data/<TF>/fnn_crossval/fold_XX/train_pos_TXX.bed

The BED has three columns (chrom, start, end), 0-based, end-exclusive.

Assumptions
----------
- FASTA headers encode genomic coordinates in the form:
      >anything::chr:start-end
  e.g.  >peak624::chrX:67537383-67537484
- Only header lines (">") are parsed; sequence lines are ignored.
- Coordinates in the header are integers.
"""
import sys
from pathlib import Path

USAGE = """\
Usage: make_train_bed.py /path/to/data

This will look under data/<TF>/fnn_crossval/fold_XX/train_pos_TXX.fa
and produce train_pos_TXX.bed alongside each train_pos_TXX.fa.
"""

if len(sys.argv) != 2:
    print(USAGE, file=sys.stderr)
    sys.exit(1)

data_root = Path(sys.argv[1])
if not data_root.is_dir():
    print(f"ERROR: “{data_root}” is not a directory", file=sys.stderr)
    sys.exit(1)

# Loop over each TF folder under data_root
for tfdir in sorted(data_root.iterdir()):
    fnn_crossval = tfdir / "fnn_crossval"
    if not fnn_crossval.is_dir():
        continue

    print(f"==== Processing TF: {tfdir.name} ====")
    # Now reuse the same logic as make_train_bed.py:
    for fold_dir in sorted(fnn_crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue

        # find train_pos_TXX.fa in this fold
        for fa_path in fold_dir.glob("train_pos_T[0-9][0-9].fa"):
            stem = fa_path.stem            # e.g. "train_pos_T01"
            bed_path = fold_dir / f"{stem}.bed"

            with open(fa_path, "r") as fin, open(bed_path, "w") as fout:
                for line in fin:
                    if not line.startswith(">"):
                        continue
                    hdr = line[1:].strip()
                    # “peak624::chrX:67537383-67537484”
                    parts = hdr.split("::", 1)
                    if len(parts) != 2 or ":" not in parts[1] or "-" not in parts[1]:
                        print(f"  WARNING: couldn’t parse {hdr}", file=sys.stderr)
                        continue
                    region = parts[1]                 # "chrX:67537383-67537484"
                    chrom, rest = region.split(":", 1)
                    start_str, end_str = rest.split("-", 1)
                    try:
                        start = int(start_str)
                        end   = int(end_str)
                    except ValueError:
                        print(f"  WARNING: non‐integer coords in {hdr}", file=sys.stderr)
                        continue
                    fout.write(f"{chrom}\t{start}\t{end}\n")
            print(f"  Wrote {bed_path}")

print("All TFs processed.")
