#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Prune TF datasets with too few peaks
#
# What this script does:
#   - For each TF under $HOME/msc_project/data/<TF>,
#     check sorted/peaks_sorted.bed (combined narrowPeak list).
#   - If it has fewer than 1800 peaks (lines), *delete*:
#       1) the entire TF data directory ($PROJ_DIR/data/<TF>)
#       2) the corresponding results directory ($PROJ_DIR/results/<TF>)
#
#    WARNING: This is destructive. It uses `rm -rf` on matched directories.
#    Consider adding a dry-run mode or a confirmation prompt before removal.
#
# Requirements:
#   - coreutils: wc, rm, basename
#   - directory layout produced by your earlier scripts
# -----------------------------------------------------------------------------

PROJ_DIR="$HOME/msc_project"                 # my base directory for ENCODE data for each TF

for tfdir in "$PROJ_DIR"/data/*; do          # loop over each TF directory in $PROJ_DIR/data/
  name=$(basename "$tfdir")                  # name = the TF directory base name, e.g. CTCF
  sorted="$tfdir/sorted/peaks_sorted.bed"    # path to the sorted peaks file for this TF
    if [ -f "$sorted" ]; then                # check if peaks_sorted.bed exists
        num_peaks=$(wc -l < "$sorted")       # count the number of lines in peaks_sorted.bed
        if [ "$num_peaks" -lt 1800 ]; then   # if the number of peaks is less than 1800
        echo "$name has $num_peaks peaks ... removing from data $tfdir and removing from results $PROJ_DIR/results/$name"
        rm -rf "$tfdir"                      # remove the TF directory
        rm -rf "$PROJ_DIR/results/$name"     # remove the results directory for this TF
        else
        fi
    fi
done
