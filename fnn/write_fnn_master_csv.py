#!/usr/bin/env python3
"""
Aggregate per-TF FNN AUC results into a single master CSV.

What this script does
---------------------
- Scans ~/msc_project/data for TF directories (one level deep).
- For each TF, looks for:
    ~/msc_project/results/<TF>/fnn/metrics/<TF>_auc.csv
  which should contain per-fold AUC values (and an 'average' row).
- Concatenates all available per-TF CSVs into one DataFrame.
- Writes the combined table to:
    ~/msc_project/master_results/fnn/metrics/master_auc.csv

Notes
-----
- The current implementation writes/updates the master CSV **inside** the loop,
  i.e., after each TF is processed. If you prefer to write it once at the end,
  move the writing block outside the `for` loop.
"""
import os
from pathlib import Path
import pandas as pd

DATA_ROOT = Path("~/msc_project/data").expanduser()
RESULT_ROOT = Path("~/msc_project/results").expanduser()
MASTER_METRICS = Path("~/msc_project/master_results/fnn").expanduser() / "metrics"
all_tf_dirs = [d for d in sorted(DATA_ROOT.iterdir()) if d.is_dir()]

master_dfs = []

for tf_dir in all_tf_dirs:
    tf_name = tf_dir.name
    auc_file = RESULT_ROOT / tf_name / "fnn/metrics" / f"{tf_name}_auc.csv"
    
    if auc_file.exists():
        try:
            df = pd.read_csv(auc_file)
            master_dfs.append(df)
            print(f"Added results for {tf_name} to master collection")
        except Exception as e:
            print(f"Warning: Could not read results for {tf_name}: {e}")
    else:
        print(f"Warning: No AUC results found for {tf_name}")

    if master_dfs:
        master_df = pd.concat(master_dfs, ignore_index=True)
        master_df.to_csv(MASTER_METRICS/f"master_auc.csv", index=False)
        print(f"Master CSV written with results from {len(master_dfs)} TFs")
    else:
        print("No results to write to master CSV")