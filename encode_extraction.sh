#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Bulk download of ENCODE optimal IDR thresholded TF ChIP-seq narrowPeak files (hg19).
#
# What this script does:
#   1) Queries the ENCODE API (JSON) for *released* hg19 narrowPeak files with
#      output_type="optimal IDR thresholded peaks" where the target is a TF.
#   2) Extracts pairs of (TF name, download URL) via jq.
#   3) Creates per-TF folders under $HOME/msc_project/data/<TF>/raw.
#   4) Downloads all peak files with wget, up to PARALLEL jobs at a time.
#
# Requirements:
#   - curl, jq, wget, bash ≥ 4 (for `wait -n`)
#   - A POSIX-like environment with job control
#
# Notes:
#   - TF names are sanitized to filesystem-safe tokens (spaces, slashes → "_").
#   - If cloud storage URLs are present, those are preferred; otherwise falls
#     back to ENCODE-hosted path (.href).
#   - Re-running continues partial downloads (`wget -c`).
# -----------------------------------------------------------------------------

set -euo pipefail    # Exit on error, undefined variable, or pipe failure

# Directory setup
BASE_DIR="$HOME/msc_project/data"  # Base directory for ENCODE data for each TF       
RAW_SUB="raw"                      # Subdirectory for raw ENCODE files (ENCFF...)                  
PARALLEL=4                         # Number of parallel downloads

# ENCODE API
ENCODE_URL='https://www.encodeproject.org/search/?type=File&assembly=hg19&file_format=bed&file_format_type=narrowPeak&output_type=optimal%20IDR%20thresholded%20peaks&target.investigated_as=transcription%20factor&status=released&limit=all&format=json'

# Create base directory if it doesn't exist 
mkdir -p "$BASE_DIR"

echo "Querying ENCODE …"

# Fetches JSON data from ENCODE API
json=$(curl -sL -H 'Accept: application/json' "$ENCODE_URL")  # -sL: silent, follow redirects, -H: set header


hits=$(jq '.total' <<<"$json") # jq to extract how many files matched the query
[[ $hits -eq 0 ]] && { echo "No files matched — exiting."; exit 1; }  # Exit if no files found
echo "Found $hits peak files. Starting downloads …"
echo

# jq emits:  <TF><TAB><URL>
while IFS=$'\t' read -r tf url; do   # read each line with tab as delimiter
    # Sanitise the TF name 
    tf=$(echo "$tf" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')           # removes leading/trailing spaces 
    safetf=$(echo "$tf" | tr '[:space:]/' '_' | tr -cs 'A-Za-z0-9_.-_' '_') # replaces unsafe characters with '_'
    safetf=${safetf##_}; safetf=${safetf%%_}                                # strip trailing '_'

    # build output path  <BASE>/<TF>/raw/ 
    dir="$BASE_DIR/$safetf/$RAW_SUB"   # raw subdirectory for each TF
    mkdir -p "$dir"                    # create directory if it doesn't exist

    echo "  $safetf <--- ${url##*/}"   # print what file is being downloaded
    (
        wget -q --show-progress -c -P "$dir" "$url"   # -q: quiet, --show-progress: show progress bar, -c: continue incomplete downloads, -P: specify directory
    ) &
    # limit number of parallel wget jobs
    (( $(jobs -r | wc -l) >= PARALLEL )) && wait -n
done < <(
    jq -r '.["@graph"][] |
           [.target.label,
            (.cloud_metadata.url // ("https://www.encodeproject.org" + .href))] |
           @tsv' <<<"$json"
)

wait   # wait for all background jobs to finish
echo -e "\n  All done. Peak files are under  $BASE_DIR/<TF>/$RAW_SUB/"