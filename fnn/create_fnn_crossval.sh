#!/bin/bash
# -----------------------------------------------------------------------------
# Prepare FNN-friendly cross-validation inputs per TF.
#
# For each transcription factor directory under ~/msc_project/data:
#   - If a crossval/ folder exists, create a parallel fnn_crossval/ folder.
#   - For folds 01..10, copy/rename the standard files:
#       crossval/fold_XX/train_TXX.fa  -> fnn_crossval/fold_XX/train_pos_TXX.fa
#       crossval/fold_XX/bg_BXX.fa     -> fnn_crossval/fold_XX/test_bg_BXX.fa
#       crossval/fold_XX/test_FXX.fa   -> fnn_crossval/fold_XX/test_pos_FXX.fa
#
# Notes:
#   - Missing files are silently ignored (cp ... 2>/dev/null).
#   - Adjust BASE_DIR if your project lives elsewhere.
#   - Add `-v` to cp for verbose copies, or remove `2>/dev/null` to see errors.
# -----------------------------------------------------------------------------

# Base data directory
BASE_DIR=~/msc_project/data

# Loop through all TF directories (excluding non-directories like files)
for TF_DIR in "$BASE_DIR"/*/; do
    # Check if crossval folder exists within this TF_DIR
    if [ -d "$TF_DIR/crossval" ]; then
        echo "Processing $(basename "$TF_DIR")..."

        # Define source and destination
        SRC_DIR="$TF_DIR/crossval"
        DEST_DIR="$TF_DIR/fnn_crossval"

        # Create destination root
        mkdir -p "$DEST_DIR"

        # Loop through folds 01 to 10
        for i in $(seq -w 1 10); do
            fold="fold_$i"
            src_fold="$SRC_DIR/$fold"
            dest_fold="$DEST_DIR/$fold"

            # Create destination fold directory
            mkdir -p "$dest_fold"

            # Copy only the desired files, if they exist
            cp "$src_fold/train_T$i.fa" "$dest_fold/train_pos_T$i.fa" 2>/dev/null
            cp "$src_fold/bg_B$i.fa" "$dest_fold/test_bg_B$i.fa" 2>/dev/null
            cp "$src_fold/test_F$i.fa" "$dest_fold/test_pos_F$i.fa" 2>/dev/null
        done
    fi
done
echo "All TFs processed."
