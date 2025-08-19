#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Concatenate and sort ENCODE TF narrowPeak files per TF by signal strength.
#
# For each TF directory under $HOME/msc_project/data:
#   - find its raw/ subfolder containing one or more *.bed / *.bed.gz peak files
#   - concatenate (and decompress if needed)
#   - sort all peaks by the narrowPeak "signalValue" (column 7, 1-based) descending
#   - write a single sorted BED: $HOME/msc_project/data/<TF>/sorted/peaks_sorted.bed
#
# Requirements:
#   - bash, find, sort, zcat (from gzip), coreutils
#   - Input files are ENCODE narrowPeak format (BED6 + 4 columns; signalValue is col 7)
#
# Notes:
#   - Keeps fields as-is; only ordering changes.
#   - Handles both .bed and .bed.gz via zcat -f.
#   - If a TF has no BEDs in raw/, that TF is skipped.
# -----------------------------------------------------------------------------

set -euo pipefail   # Exit on error, undefined variable, or pipe failure

# Directory setup 
BASE_DIR="$HOME/msc_project/data"   # my base directory for ENCODE data for each TF
RAW_SUB="raw"                       # BASE_DIR/<TF>/raw/ will contain the downloaded peak files
OUT_SUB="sorted"                    # output directory for peaks_sorted.bed.gz here
SIGNAL_COL=7                        # narrowPeak file: signalValue column (1‑based)

echo "Scanning $BASE_DIR for TF peak folders …"

# Walk every <BASE_DIR>/<TF>/raw/
find "$BASE_DIR" -type d -name "$RAW_SUB" | while read -r RAWDIR; do
    TFDIR=$(dirname "$RAWDIR")          # gets the parent directory of <raw> e.g. $BASE_DIR/CTCF
    OUTDIR="$TFDIR/$OUT_SUB"            # output directory for sorted peaks
    mkdir -p "$OUTDIR"                  # create output directory if it doesn't exist

    # List BED peak files in this RAWDIR
    mapfile -t FILES < <(ls "$RAWDIR"/*.bed "$RAWDIR"/*.bed.gz 2>/dev/null || true)
    # Skip if no BED files found
    if (( ${#FILES[@]} == 0 )); then
        echo "No BEDs found in $RAWDIR — skipping."
        continue
    fi

    TMP="$OUTDIR/tmp.concat.bed"         # path for temporary concatenated file

    if (( ${#FILES[@]} == 1 )); then
        echo "--> Single replicate in $RAWDIR"
        zcat -f "${FILES[0]}" > "$TMP"   # if only one file, just decompress and copy to $TMP
    else
        echo "→ Concatenating ${#FILES[@]} replicates in $RAWDIR"
        zcat -f "${FILES[@]}" > "$TMP"   # if multiple files, concatenate, decompress, and copy to TMP
    fi

    # Sort by signalValue (col 7) descending
    sort -k${SIGNAL_COL},${SIGNAL_COL}nr "$TMP" \   # sort concatenated BED file by signalValue, numeric, descending order
      > "$OUTDIR/peaks_sorted.bed"                  # write to peaks_sorted.bed file   

    rm "$TMP"   # remove temporary file
    echo "->  Wrote  $OUTDIR/peaks_sorted.bed"
done

echo "All datasets processed."