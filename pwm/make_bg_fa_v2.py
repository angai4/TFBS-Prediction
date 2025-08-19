#!/usr/bin/env python3
"""
GC-matched background FASTA generator with fallback.

What this script does:
  1) Reads a precomputed background window file (TSV) containing:
       chrom, start, end, GC%, sequence
     and groups windows into integer GC% bins (0–100).
  2) For each TF's cross-validation fold, reads the test-set GC summary
     (`bg_BXX.gc.txt`) to count how many background windows are needed per GC% bin.
  3) Samples background windows per GC bin to match the test-set distribution.
     If a bin is short on windows, falls back to lower GC bins (target-1, -2, ...)
     up to a configurable distance (default: 50%).
  4) Writes a FASTA file with headers `>chr:start-end` and the sampled sequences.

Assumptions / inputs:
  - Project layout:
      ~/msc_project/data/<TF>/crossval/fold_XX/bg_BXX.gc.txt   (per-fold GC summary; chrom start end GC%)
      ~/msc_project/background/windows.gc.seq.txt              (global background: chrom start end GC% seq)
  - The background file has whitespace-delimited columns (tab or spaces).
  - GC values in both files are interpreted as integer percentages (0–100).

Outputs:
  - Per fold: `~/msc_project/data/<TF>/crossval/fold_XX/bg_BXX.fa` (FASTA with matched-GC background sequences)

Notes / caveats:
  - Sampling is random and not seeded → results vary run-to-run; set `random.seed(...)` if you need reproducibility.
  - Fallback only searches *downward* in GC (target-1, target-2, ...). If you want symmetric fallback,
    extend the strategy to try ± offsets.
  - The script does not prevent reusing the same background window across *different* bins/folds;
    if you require global uniqueness, track and exclude already-selected coordinates.
"""

import sys, random
from pathlib import Path
from collections import Counter, defaultdict

DATA_ROOT = Path.home() / "msc_project" / "data" # Path to the data directory
BG_ROOT = Path.home() / "msc_project" / "background" / "windows.gc.seq.txt" # Path to the background windows file

# Group background windows by GC-bin
print(f"Collecting 'chr, start, end, seq' for each GC% in {BG_ROOT} ...")
bg_by_bin = defaultdict(list)   # dict[int GC%] --> list of (chr, start, end, seq)
with open(BG_ROOT) as g:
    for line in g:
        line = line.strip()
        if not line or line.startswith('#'):
            # skip empty lines and comments
            continue
        # expect at least 5 fields: chrom, start, end, gc, seq
        # split into 5 pieces max so the sequence remains intact
        chr_, start, end, gc, seq = line.split(None, 4)
        # bin GC by integer percentage; "70.0" becomes 70
        bin_gc = int(float(gc))
        bg_by_bin[bin_gc].append((chr_, start, end, seq))
print(f"Collected {len(bg_by_bin.keys())} GC% bins from background data")

def find_fallback_sequences(target_gc, needed_count, bg_by_bin, max_fallback_distance=50):
    """
    Find sequences for a given GC%, falling back to n-1 GC% if insufficient data.
    
    Args:
        target_gc: Target GC percentage
        needed_count: Number of sequences needed
        bg_by_bin: Dictionary of sequences by GC bin
        max_fallback_distance: Maximum distance to fall back (default 50%)
    
    Returns:
        List of selected sequences, or None if insufficient data found
    """
    for gc_offset in range(max_fallback_distance + 1):
        # Try target GC% first, then target-1, target-2, etc.
        fallback_gc = target_gc - gc_offset
        
        # Skip negative GC values
        if fallback_gc < 0:
            continue
            
        pool = bg_by_bin.get(fallback_gc, [])
        
        # if the pool has enough windows, sample without replacement (within this bin)
        if len(pool) >= needed_count:
            if gc_offset > 0:
                print(f"  Fallback: Using {fallback_gc}% GC (needed {target_gc}% but insufficient data)")
            return random.sample(pool, needed_count)
    
    return None  # Insufficient data even with fallback

# Loop over each TF folder under data_root
for tfdir in sorted(DATA_ROOT.iterdir()):
    crossval = tfdir / "crossval"
    if not crossval.is_dir():
        continue
        
    print(f"==== Processing TF: {tfdir.name} ====\n")
    # iterate over each fold directory
    for fold_dir in sorted(crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue
        
        print(f"Processing: {fold_dir.name}")

        # Count how many test windows fall into each GC% bin for this fold
        test_counts = Counter()
        print(f"Resetting GC content counts: \n{test_counts}")
        fold_number = fold_dir.name.split('_')[-1]
        name = f"bg_B{fold_number.zfill(2)}"
        
        # the per-fold GC summary file created by the previous script
        with open(fold_dir / f"{name}.gc.txt") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # expect four fields: chrom, start, end, gc
                chr_, start, end, gc = line.split()
                bin_gc = int(float(gc))
                test_counts[bin_gc] += 1
        
        print(f"{tfdir.name} counting GC content in {name}.gc.txt: \n{len(test_counts)} GC% bins found")

        # Sample with fallback strategy
        selected = []
        print(f"Resetting the list of selected sequences ----> {selected}") # 0
        skip_tf = False # Flag to skip TF if insufficient background data even with fallback
        
        for bin_gc, want in test_counts.items():
            fallback_sequences = find_fallback_sequences(bin_gc, want, bg_by_bin)
            
            if fallback_sequences is None:
                print(f"ERROR: need {want} windows at {bin_gc}% GC but insufficient data even with fallback")
                print("Skipping this TF due to insufficient background data (even with fallback)")
                skip_tf = True
                break
            
            selected.extend(fallback_sequences)
            
        if skip_tf:
            break # Skip to next TF if insufficient data even with fallback
            
        print(f"{len(selected)}: sequences in the selected sequences list from background data for {fold_dir} at GC% bins") # Should be the same as the amount of test sequences
        
        # Write out FASTA with headers >chr:start-end
        with open(fold_dir / f"{name}.fa", 'w') as out:
            for chr_, start, end, seq in selected:
                out.write(f">{chr_}:{start}-{end}\n{seq}\n")
        print(f"Wrote {len(selected)} bg sequences to {fold_dir / f'{name}.fa'}")
    print(f"Completed processing for TF: {tfdir.name}\n")

print("All TFs processed successfully.")