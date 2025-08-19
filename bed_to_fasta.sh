#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Select top-N peaks and extract 101 bp summit-centered sequences per TF.
#
# For each TF with a sorted peak file at:
#   $HOME/msc_project/data/<TF>/sorted/peaks_sorted.bed
# this script:
#   1) Takes the top 600 rows (already sorted by signalValue descending).
#   2) Converts each narrowPeak row to a 101 bp window centered on the summit:
#        - If column 10 (summit offset) == -1, use the region midpoint.
#        - Else, summit = start (col 2) + summit_offset (col 10).
#        - Window = [summit-50, summit+50], written as BED (0-based, end-exclusive).
#   3) Extracts the corresponding sequences from hg19 into FASTA
#      for motif discovery (e.g., MEME).
#
# Requirements:
#   - bedtools (getfasta)
#   - samtools (to create/ensure <genome>.fai exists)
#   - Input peaks_sorted.bed must be narrowPeak format (BED6 + 4 extra cols).
#
# Notes:
#   - The genome FASTA must be indexed: a sibling file <hg19.fa>.fai is required.
#   - This script ignores strand when extracting sequences (suitable for MEME).
#   - peaks_sorted.bed is assumed to be sorted by signalValue (col 7) descending.
# -----------------------------------------------------------------------------

set -euo pipefail   # exit on error, undefined variable, or pipe failure
shopt -s nullglob   # unmatched globs will expand to an empty array

# Directory setup
BASE_DIR="$HOME/msc_project/data"   # my base directory for ENCODE data for each TF
GENOME="$HOME/genomes/hg19.fa"      # reference FASTA (indexed with samtools faidx)
TOP=600                             # number of strongest peaks to keep
WIN=50                              # bp flanking each side of peak summit (=> 101‑bp window)
OUT_SUB="top600"                    # output folder name

# Check for genome index
[[ -f ${GENOME}.fai ]] || { echo "${GENOME}.fai index not found"; exit 1; }

echo "Creating ±${WIN} bp windows around peak‑max for top $TOP peaks …"    # print message to user

# loop over every peaks_sorted.bed
find "$BASE_DIR" -name peaks_sorted.bed -path "*/sorted/*" -print0 |   # searches for files named peaks_sorted.bed specifically in */sorted/* directories
while IFS= read -r -d '' BEDFILE; do                                   # BEDFILE is the full path to peaks_sorted.bed
    TFDIR=$(dirname "$(dirname "$BEDFILE")")                           # get the parent directory of peaks_sorted.bed, e.g. $BASE_DIR/<TF>/sorted
    OUTDIR="$TFDIR/$OUT_SUB"                                           # e.g. $BASE_DIR/<TF>/top600
    mkdir -p "$OUTDIR"                                                 # create output directory if it doesn't exist

    echo "$(basename "$TFDIR")"   # print TF name

    # 1) extract top $TOP (=600) rows and output to $BASE_DIR/<TF>/top600/top600.narrowPeak
    head -n $TOP "$BEDFILE" > "$OUTDIR/top${TOP}.narrowPeak"         # BEDFILE is the full path to peaks_sorted.bed

    # 2) convert each row to a 101bp centred window
    awk -v w=$WIN 'BEGIN{OFS="\t"}                                   # -v means to pass a variable to awk, the command line variable w (50) is passed to awk as a variable with the same name
        {
          summit = ($10 == -1 ? int(($2+$3)/2) : $2 + $10);          # if summit offset is -1, use midpoint of $2 and $3, otherwise use $2 + $10
          start  = summit - w;                                       # start of the window is summit - w (50)
          if (start < 0) start = 0;                                  # BED start cannot be negative, so set it to 0 if it is
          end    = summit + w + 1;                                   # end of the window is summit + w + 1 (to include the base at end)
          print $1, start, end, "peak"NR;                            # print chromosome, start, end, and a name for the peak
        }' "$OUTDIR/top${TOP}.narrowPeak" \                          # input file is the top600.narrowPeak file
        > "$OUTDIR/top${TOP}_$((WIN*2+1))bp.bed"                     # output file is top600_101bp.bed

    # 3) FASTA extraction with bedtools getfasta (ignore strand for MEME input)
    bedtools getfasta -fi "$GENOME" \                                 # -fi   <input FASTA file>
                      -bed "$OUTDIR/top${TOP}_$((WIN*2+1))bp.bed" \   # -bed  <input BED file>
                      -fo  "$OUTDIR/top600_fasta.fa" \                # -fo   <output file name>
                      -name                                           # -name  use the name field and coordinates for the FASTA header

    echo "->  $OUTDIR/top600_fasta.fa  ($(grep -c '^>' "$OUTDIR/top600_fasta.fa") seqs)"   # print number of sequences in the FASTA file
done

echo "Finished – FASTAs ready for MEME in each  $OUT_SUB/  folder."   # print final message to user