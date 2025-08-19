#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Run MAST motif scanning on per-fold training FASTAs for each TF.
#
# For each TF under $HOME/msc_project/data:
#   - locate its cross-validation folds (…/crossval/fold_XX/)
#   - index the fold’s training FASTA (train_TXX.fa) with samtools faidx
#   - run MAST using the TF-specific MEME motif (results/<TF>/motif/meme.meme)
#   - parse MAST hits into BED, then extract hit sequences as FASTA
#
# Requirements:
#   - MEME Suite installed and `mast` available on PATH
#   - samtools available on PATH (for `samtools faidx`)
#   - bedtools available on PATH (for `bedtools getfasta`)
#   - Directory layout:
#       $HOME/msc_project/data/<TF>/crossval/fold_XX/train_TXX.fa
#       $HOME/msc_project/results/<TF>/motif/meme.meme   (from a prior MEME run)
#
# Notes:
#   - The AWK parser assumes a specific MAST tabular output layout:
#       $1=seq_id, $2=orientation (+1/-1), $5=start(1-based), $6=end(1-based)
#     Adjust field numbers if your MAST version outputs a different column order.
#   - If a fold produces ≤3 lines in mast.txt (typically only headers), the script
#     skips the *entire* TF (using `continue 2`).
# -----------------------------------------------------------------------------
set -euo pipefail   # exit on error, undefined variable, or pipe failure
shopt -s nullglob   # unmatched globs will expand to an empty array

# Directory setup
PROJ_DATA="$HOME/msc_project/data"            # my base directory for ENCODE data for each TF    
PROJ_RESULTS="$HOME/msc_project/results"      # results directory
MT="1e-4"                                     # p‐value threshold
EV="10"                                       # e‐value threshold


for TFDIR in "$PROJ_DATA"/*; do               # loop over each TF directory in $PROJ_DIR/data/
  [ -d "$TFDIR" ] || continue                 # skip anything thats not a directory
  TF=$(basename "$TFDIR")                     # TF = the TF directory base name, e.g. CTCF
  XCVAL="$TFDIR/crossval"                     # crossval directory
  MOTIF="$PROJ_RESULTS/$TF/motif/meme.meme"   # where the meme output lives

  # skips the TF if its motif file or cross-validation directory is missing
  if [ ! -f "$MOTIF" ]; then
    echo "  Motif for $TF not found: $MOTIF" >&2
    continue
  fi
  if [ ! -d "$XCVAL" ]; then
    echo "  crossval dir for $TF not found: $XCVAL" >&2
    continue
  fi

  echo " Processing TF $TF"
  for FOLD in "$XCVAL"/fold_*; do       # iterates over each cross-validation fold directory
    [ -d "$FOLD" ] || continue          # skips anything thats not a directory
    FOLDNAME=$(basename "$FOLD")        # e.g. FOLDNAME = fold_01
    NUM=${FOLDNAME#fold_}               # e.g. NUM = 01

    echo "  ->  $FOLDNAME"
    TRAIN="$FOLD/train_T${NUM}.fa"
    if [ ! -f "$TRAIN" ]; then
      echo "      Missing training FASTA: $TRAIN" >&2
      continue
    fi

    OUT="$FOLD/mast_out"
    mkdir -p "$OUT"

    # 1) index the FASTA
    samtools faidx "$TRAIN"

    # 2) run MAST
    echo "     Running MAST on $TRAIN …"
    mast \
      -oc       "$OUT" \
      -hit_list \
      -best     \
      -nostatus \
      -mt "$MT" \
      -ev "$EV" \
      "$MOTIF" \
      "$TRAIN" \
      > "$OUT/mast.txt"

    # 2b) if mast.txt has ≤3 lines (header only), skip entire TF
    lines=$(wc -l < "$OUT/mast.txt")
    if [ "$lines" -le 3 ]; then
      echo "      No real MAST hits (only $lines header lines) → skipping TF $TF"
      continue 2
    fi

    # 3) parse mast.txt → mast_hits.bed
    awk 'BEGIN { OFS="\t" }
         /^#/   { next }
         {
           id          = $1
           start0      = $5 - 1
           if (start0 < 0) start0 = 0
           end1        = $6
           strand_sign = ($2 == "+1" ? "+" : "-")
           strand_num  = 1
           name        = id ":" start0 "-" end1
           print id, start0, end1, name, strand_num, strand_sign
         }' \
      "$OUT/mast.txt" > "$OUT/mast_hits.bed"

    # 4) extract subsequences
    bedtools getfasta \
      -fi  "$TRAIN" \
      -bed "$OUT/mast_hits.bed" \
      -nameOnly -s \
      -fo  "$OUT/mast_hits.fa"

    COUNT=$(grep -c '^>' "$OUT/mast_hits.fa")
    echo "     extracted $COUNT hits → mast_hits.fa"
  done

  echo " Completed TF $TF"
done

echo " All TFs processed!"
