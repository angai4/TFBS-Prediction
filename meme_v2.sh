#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# MEME motif discovery runner for ENCODE-style TF folders.
#
# For each TF directory under $HOME/msc_project/data, this script:
#   1) finds the FASTA file at .../<TF>/top600/top600_fasta.fa
#   2) runs MEME with DNA mode, ZOOPS, one motif (width 6–20), using N cores
#   3) saves outputs under $HOME/msc_project/results/<TF>/motif/meme_out
#   4) copies MEME’s motif file (meme.txt) to results/<TF>/motif/meme.meme
#
# Requirements:
#   - MEME Suite (`meme`) on PATH
#   - Directory layout: $HOME/msc_project/data/<TF>/top600/top600_fasta.fa
#
# Notes:
#   - Re-running may reuse/overwrite existing output directories/files.
#   - Adjust MEME_CORES/MEME_ARGS as needed for your cluster/machine.
# -----------------------------------------------------------------------------
set -euo pipefail                                         # Exit on error, undefined variable, or pipe failure
shopt -s nullglob                                         # unmatched globs will expand to an empty array

# Directory setup
PROJ_DIR="$HOME/msc_project"                              # my base directory for ENCODE data for each TF
MEME_CORES=8                                              # number of cores to use for MEME
MEME_ARGS="-dna -mod zoops -nmotifs 1 -minw 6 -maxw 20"   # MEME arguments: -dna for DNA sequences, -mod zoops for zero or one occurrence per sequence, -nmotifs 1 to find one motif, -minw 6 and -maxw 20 for motif width range   

# 1) Gather all of the FASTA inputs into an array
# uses find to locate all top600_fasta.fa files under $PROJ_DIR/data
# mapfile reads the output into an array 'FASTAS', splitting on null characters
mapfile -d '' FASTAS < <(            
  find "$PROJ_DIR/data" \
       -type f \
       -name "top600_fasta.fa" \         
       -print0  
) 

# Check if any FASTA files were found
if (( ${#FASTAS[@]} == 0 )); then
  echo "No top600_fasta.fa files found under $PROJ_DIR/data"
  exit 1
fi

# Print the number of FASTA files found
echo "Found ${#FASTAS[@]} FASTA files. Beginning MEME runs…"

# 2) Loop over them in the current shell
for FASTA in "${FASTAS[@]}"; do
    TOPDIR=$(dirname "$FASTA")                     # …/<TF>/top600
    TFDIR=$(dirname "$TOPDIR")                     # …/<TF>
    TF_ID=$(basename "$TFDIR")                     # <TF>

    RESDIR="$PROJ_DIR/results/$TF_ID/motif"        # setup results directory
    OUTDIR="$RESDIR/meme_out"                      # output directory for MEME results
    mkdir -p "$OUTDIR"                             # create output directory if it doesn't exist

    # Print the current TF and FASTA file being processed
    echo "MEME on  $TF_ID  (input: $FASTA)"
    # 3) if MEME fails for one TF, catch & continue
    # run MEME with the specified arguments
    if ! meme "$FASTA" $MEME_ARGS -oc "$OUTDIR" -p "$MEME_CORES"; then
      echo "MEME failed for $TF_ID, skipping."
      continue
    fi

    cp -f "$OUTDIR/meme.txt" "$RESDIR/meme.meme"   # copy the MEME output to the results directory
    echo "Motif ready --> $RESDIR/meme.meme"       # print the output file path
done

echo "All done."                           
