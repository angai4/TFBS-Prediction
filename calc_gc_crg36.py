#!/usr/bin/env python3
# This script builds a background set of uniquely mappable genomic windows from a CRG mappability bigWig.
# Pipeline overview:
#   1) Convert the CRG 36mer mappability bigWig -> BedGraph.
#   2) Keep only bases/intervals with uniqueness == 1 (fully uniquely mappable).
#   3) Sort & merge adjacent/overlapping uniquely mappable regions.
#   4) Keep merged regions that are ≥202 bp (so they can hold at least two 101 bp windows).
#   5) Tile the remaining regions with non-overlapping 101 bp windows (step = 101).
#   6) Compute sequence features with `bedtools nuc` on hg19 (GC%, sequence).
#   7) Extract chr, start, end, GC% (0–100), and the sequence into a compact TSV.
#
# Requirements on PATH / data:
#   - UCSC `bigWigToBedGraph`
#   - `bedtools` (providing `merge`, `makewindows`, and `nuc`)
#   - `awk`, `sort`
#   - Reference FASTA: hg19.fa (and its .fai index alongside) under CRG_ROOT
#   - Input bigWig: wgEncodeCrgMapabilityAlign36mer.bigWig under CRG_ROOT
#
# Notes:
#   - The script assumes bedtools `nuc -seq` emits GC fraction and sequence in the columns referenced below.
#     Column positions may vary across versions; adjust indices if needed.

import os, sys
import subprocess
from pathlib import Path

# Root directories
PROJ_DIR = Path.home() / "msc_project"    # project root under the user's home directory
BG_ROOT = PROJ_DIR / "background"         # output directory for background windows
CRG_ROOT = Path.home() / "genomes"        # location of wgEncodeCrgMapabilityAlign36mer.bigWig and hg19.fa

# sanity checks
if not os.path.exists(BG_ROOT):
    # ensure the background directory exists
    print(f"Background directory {BG_ROOT} does not exist.\nCreating it...")
    os.makedirs(BG_ROOT, exist_ok=True)

print("Converting bigWig to BedGraph...")
if not CRG_ROOT.exists():
    # fail early if the genomes directory does not exist
    print(f"CRG root directory {CRG_ROOT} does not exist.")
    sys.exit(1)

# convert CRG 36mer mappability bigWig to BedGrapgh
# BedGraph columns: chrom, start, end, value(=uniqueness score)
subprocess.run([
    'bigWigToBedGraph',
    str(CRG_ROOT / "wgEncodeCrgMapabilityAlign36mer.bigWig"),
    str(BG_ROOT / "crg36.map.bedGraph")
])

print("===bigWig to BedGraph conversion completed===")

print("------>BedGraph Filtering<------\
\nFiltering for sequence uniqueness == 1...")
# keep only intervals whose mappability score is 1
# output a 3-column BED (chrom, start, end)
cmd1 = f"awk '$4 == 1 {{ print $1, $2, $3 }}' OFS='\t' {BG_ROOT / 'crg36.map.bedGraph'} > {BG_ROOT / 'crg36.map.bed'}"
subprocess.run(cmd1, shell=True, check=True)

print("Sorting and merging intervals in bed file...")
# sort by chrom and start, then merge overlapping intervals to create maximal uniquely mappable blocks
cmd2 = f"sort -k1,1 -k2,2n {BG_ROOT / 'crg36.map.bed'} | bedtools merge -i - > {BG_ROOT / 'crg36.map.merged.bed'}"
subprocess.run(cmd2, shell=True, check=True)

print("Filtering merged bed file for intervals >=202bp...")
# retain only merged intervals at least 202bp long so they contain at least 2 101bp windows without overlap
cmd3 = f"awk '($3 - $2) >= 202' OFS='\t' {BG_ROOT / 'crg36.map.merged.bed'} > {BG_ROOT / 'crg36.map.merged.202bp.bed'}"
subprocess.run(cmd3, shell=True, check=True)

print("Creating sliding windows of 101bp with 101bp step...")
# tile each retained interval with non-overlapping 101bp windows (step size = window size)
# the trailing awk guarantees exact 101 bp windows
cmd4 = f"bedtools makewindows -b {BG_ROOT / 'crg36.map.merged.202bp.bed'} -w 101 -s 101 | awk '$3 - $2 == 101' > {BG_ROOT / 'crg36.map.merged.101bp.windows.bed'}"
subprocess.run(cmd4, shell=True, check=True)

print("===filtering completed===")

print("Running bedtools nuc for GC% calculation...")
# compute nucleotide composition and GC% for each 101bp window against hg19
# the -seq flag appends the actual sequence to the output
cmd = f"bedtools nuc -fi {CRG_ROOT / 'hg19.fa'} -bed {BG_ROOT / 'crg36.map.merged.101bp.windows.bed'} -seq > {BG_ROOT / 'windows.nuc.txt'}"
subprocess.run(cmd, shell=True, check=True)
print("===bedtools nuc completed===")

print("Extracting chr, start, end, GC%, seq from nuc output...")
# parse bedtools nuc output, skipping comment lines. Extract:
# chrom (col 1), start (col 2), end (col 3), 
# GC% (col 5) --> convert to integer percent,
# sequence (col 13)
with open(BG_ROOT / "windows.nuc.txt") as f, open(BG_ROOT / "windows.gc.seq.txt", "w") as out:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) < 13:  # Ensure there are enough columns
            continue
        chrom, start, end, gc, seq = parts[0], parts[1], parts[2], int(round(float(parts[4]) * 100)), parts[12]
        out.write(f"{chrom}\t{start}\t{end}\t{gc}\t{seq}\n")
print("===GC% and sequence extraction completed===")

