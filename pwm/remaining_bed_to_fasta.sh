#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Hold-out complement: build 101 bp summit-centered FASTAs from all peaks
# EXCEPT the top N (default N=600) for every TF.
#
# For each TF under $HOME/msc_project/data/<TF>:
#   1) Take peaks_sorted.bed (narrowPeak, sorted by signalValue desc).
#   2) Drop the first TOP lines (keep lines TOP+1 … end).
#   3) Convert each remaining peak to a ±WIN bp window around its summit
#      (col10 = summit offset; if -1, fall back to midpoint).
#   4) Extract sequences from hg19 (FASTA) to crossval/remaining_fasta.fa.
#
# Results are used later for background estimation / cross-validation splits.
#
# Requirements:
#   - bedtools (getfasta)
#   - Reference FASTA present at $HOME/genomes/hg19.fa, indexed alongside (.fai).
#   - Input peaks are ENCODE narrowPeak format (BED6 + 4 cols; summit in col 10).
# -----------------------------------------------------------------------------

set -euo pipefail                           # exit on error, undefined variable, or pipe failure
shopt -s nullglob                           # unmatched globs will expand to an empty array

# Directory setup
PROJ_DIR="$HOME/msc_project"                # my base directory for ENCODE data for each TF
GENOME="$HOME/genomes/hg19.fa"              # reference FASTA (indexed with samtools faidx)
TOP=600                                     # number of peaks previously held out for motif discovery
WIN=50                                      # number of bps to flank on each side of summit

for tfdir in "$PROJ_DIR"/data/*; do         # loop over each TF directory in $PROJ_DIR/data/
  name=$(basename "$tfdir")                 # name = the TF directory base name, e.g. CTCF 
  echo "=== $name ==="                      # print the TF name 
  sorted="$tfdir/sorted/peaks_sorted.bed"   # path to the sorted peaks file for this TF
  cvdir="$tfdir/crossval"                   # path to the cross-validation directory for this TF
  mkdir -p "$cvdir"                         # create the cross-validation directory if it doesn't exist

  # 1) grab all but the top $TOP peaks i.e., from 601 to the end
  echo "selecting peaks #$((TOP+1)) … end"  
  cat "$sorted" \                           # read the sorted peaks file
    | tail -n +$((TOP+1)) \                 # pipe to tail to skip the first $TOP lines 
    > "$cvdir"/remaining.bed                # output to remaining.bed in the cross-validation directory

  # 2) make ±WIN windows around summit (col10 = summit offset)
  awk -v w=$WIN -v offset=$TOP 'BEGIN{OFS="\t"}         # -v w is the window size (50), -v offset is the number of peaks previously held out (600)
  {
    summit = ($10 == -1 ? int(($2+$3)/2) : $2 + $10);   # if summit offset is -1, use midpoint of $2 and $3, otherwise use $2 + $10
    start  = summit - w;                                # start of the window is summit - w (50)
    if(start<0) start=0;                                # BED start cannot be negative, so set it to 0 if it is
    end    = summit + w + 1;   # BED end is exclusive   # end of the window is summit + w + 1 (to include the base at end)
    name   = "peak" (offset + NR);                      # name for the peak, e.g. peak601, peak602, etc.
    print $1, start, end, name;                         # print chromosome, start, end, and the name for the peak
  }' "$cvdir"/remaining.bed \                           # input file is remaining.bed
    > "$cvdir"/remaining_$((WIN*2+1))bp.bed             # output file is remaining_101bp.bed

  # 3) extract FASTA of all remaining windows
  echo "extracting FASTA for $(( $(wc -l <"$cvdir"/remaining_$((WIN*2+1))bp.bed) )) peaks"
  bedtools getfasta \                                    # use bedtools getfasta to extract sequences
    -fi "$GENOME" \                                      # -fi <input FASTA file>
    -bed "$cvdir"/remaining_$((WIN*2+1))bp.bed \         # -bed <input BED file>
    -name \                                              # -name  use the name field and coordinates for the FASTA header
    -fo "$cvdir"/remaining_fasta.fa                      # -fo <output file name>

  count=$(grep -c '^>' "$cvdir"/remaining_fasta.fa)      # count the number of sequences in the FASTA file
  echo "-->  $cvdir/remaining_fasta.fa  ($count seqs)"   # print the output file path and number of sequences
done

echo "Finished – FASTAs ready for splitting into train/test in each crossval/ folder."