#!/usr/bin/env python3
import sys
from pathlib import Path

USAGE = """\
Usage: make_test_bed_all.py /path/to/data

This will look under data/<TF>/crossval/fold_XX/test_FXX.fa
and produce test_FXX.bed alongside each test_FXX.fa.
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
    crossval = tfdir / "crossval"
    if not crossval.is_dir():
        continue

    print(f"==== Processing TF: {tfdir.name} ====")
    # Now reuse the same logic as make_test_bed.py:
    for fold_dir in sorted(crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue

        # find test_FXX.fa in this fold
        for fa_path in fold_dir.glob("test_F[0-9][0-9].fa"):
            stem = fa_path.stem            # e.g. "test_F01"
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
