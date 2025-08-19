#!/usr/bin/env python3
"""
Compute PWM-based ROC/AUC per TF and per cross-validation fold, then save plots and CSV metrics.

Pipeline (for each TF under ~/msc_project/data/<TF>):
  1) Estimate background nucleotide frequencies (A/C/G/T) from crossval/remaining_fasta.fa.
  2) For each fold (crossval/fold_XX):
       - Read aligned motif instances from mast_out/mast_hits.fa.
       - Build a Biopython motif from those sequences and derive a log-odds PWM
         using the empirically estimated background frequencies (with pseudocounts).
       - Score both the positive set (test_FXX.fa) and the GC-matched background (bg_BXX.fa)
         with the PWM and its reverse complement; for each sequence keep the max score.
       - Compute ROC and AUC (scikit-learn), plot ROC to PNG, and record metrics.
  3) Write per-TF AUCs (per fold + an average row) to results/<TF>/pwm/metrics/<TF>_auc.csv.
  4) Aggregate all TFs into master_results/pwm/metrics/master_auc.csv (overwritten each TF iteration).

Requirements:
  - Biopython (SeqIO, motifs)
  - scikit-learn (roc_curve, roc_auc_score)
  - matplotlib
  - Directory layout produced by your prior steps (MAST outputs, test/background FASTAs).
"""

import os
from pathlib import Path
from Bio import SeqIO, motifs
from collections import Counter
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
import csv

# Root directories
PROJECT = Path.home() / "msc_project"
DATA_ROOT = PROJECT / "data"
RESULTS_ROOT = PROJECT / "results"

# Track all AUCs for the master file
master_auc_records = []

# Loop over each TF folder under data_root
for tfdir in sorted(DATA_ROOT.iterdir()):
    crossval = tfdir / "crossval"
    if not crossval.is_dir():
        # skip TFs lacking the expected crossval layout
        continue
        
    print(f"\n====== Processing: {tfdir.name} ======")

    # Initialise AUC records for this TF
    auc_records = []

    # Calculate actual background frequencies from the remaining_fasta.fa file
    rem_fasta = crossval / "remaining_fasta.fa"
    if not rem_fasta.is_file():
        print(f"\nSkipping  {tfdir.name}: No remaining_fasta.fa file found.")
        continue

    counts = Counter()
    # count A/C/G/T across all sequences (case-insensitive), ignoring any other characters
    for rec in SeqIO.parse(rem_fasta, "fasta"):
        counts.update(b for b in str(rec.seq).upper() if b in "ACGT")

    # convert counts to frequencies
    total = sum(counts[b] for b in "ACGT")
    background = {b: counts[b]/total for b in "ACGT"}
    print(f"\nBackground freqeuncies: {background}")
    
    # Loop over each fold directory
    for fold_dir in sorted(crossval.glob("fold_[0-9][0-9]")):
        if not fold_dir.is_dir():
            continue
        
        print(f"\nProcessing: {fold_dir.name}")

        # Look for aligned MAST sequences
        hits_fa = fold_dir / "mast_out" / "mast_hits.fa"
        if not hits_fa.is_file():
            print(f"  Skipping fold {fold_dir.name}: No mast_hits.fa file found.")
            continue

        # Read in aligned MAST sequences
        sequences = [str(rec.seq).upper() for rec in SeqIO.parse(hits_fa, "fasta")]
        print("  Aligned sequences:", len(sequences))

        # Build a motif object
        m = motifs.create(sequences)
        print("  Motif length:", len(m))

        # Create a PWM (bg is 0.25 by default), we use the actual background frequencies
        pwm = m.counts.normalize(pseudocounts=1).log_odds(background)
        pwm_rc = pwm.reverse_complement()
        print("  PWM and reverse complement PWM created.")
        
        # Function to get max score (fwd vs rev) for each sequence
        def get_max_scores(fasta_path):
            """
            Slide the PWM across each sequence and return the max log-odds score,
            considering both forward and reverse-complement PWMs.
            Skips sequences shorter than the motif or containing non-ACGT bases.
            """
            scores = []
            for rec in SeqIO.parse(fasta_path, "fasta"):
                seq = str(rec.seq).upper()

                # discard if shorter than motif width
                if len(seq) < len(pwm):
                    continue

                # discard ambiguous sequences
                if any(base not in "ACGT" for base in seq):
                    print(f"  Skipping {rec.id}: contains non-ACGT bases")
                    continue

                try:
                    fwd_score = max(pwm.calculate(seq))
                    rev_score = max(pwm_rc.calculate(seq))
                    score = max(fwd_score, rev_score)

                    if not (score != score): # Check for NaN (Nan != Nan is True)
                        scores.append(score)
                except Exception as e:
                    print(f"  Skipping {rec.id} due to scoring error: {e}")
                    continue
            return scores
        
        # Get test and background fasta files
        fold_number = fold_dir.name.split('_')[-1]
        test_name = f"test_F{fold_number.zfill(2)}.fa"
        bg_name = f"bg_B{fold_number.zfill(2)}.fa"

        test_seqs = fold_dir / test_name
        bg_seqs = fold_dir / bg_name
    
        # Get PWM scores
        print("  Calculating max PWM scores for test and background sequences...")
        pos_scores = get_max_scores(test_seqs)
        neg_scores = get_max_scores(bg_seqs)
        y_scores = pos_scores + neg_scores
        y_true = [1] * len(pos_scores) + [0] * len(neg_scores)

        # Calculate ROC and AUC
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        auc = roc_auc_score(y_true, y_scores)
        print(f"  AUC: {auc:.3f}")

        # record per-fold AUC for this TF
        auc_records.append({
            "TF": tfdir.name,
            "Fold": fold_dir.name,
            "AUC": round(auc, 3)
        })

        # Plot and save ROC
        fig_dir = RESULTS_ROOT / tfdir.name / "pwm" / "figures" 
        fig_dir.mkdir(parents=True, exist_ok=True)
        roc_path = fig_dir / f"{tfdir.name}_{fold_dir.name}_roc.png"
        
        plt.figure()
        plt.plot(fpr, tpr, label=f"ROC curve (AUC = {auc:.3f})")
        plt.plot([0, 1], [0, 1], 'k--', label="Random Classifier")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"{tfdir.name} - {fold_dir.name} PWM ROC")
        plt.legend(loc="lower right")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{roc_path}", dpi=300)
        plt.close()
        print(f"  ROC curve saved as: {roc_path}")

    # Calculate and append average AUC
    if auc_records:
        mean_auc = round(sum(r["AUC"] for r in auc_records) / len(auc_records), 3)
        auc_records.append({
            "TF": tfdir.name,
            "Fold": "Average",
            "AUC": mean_auc
        })
        print(f"\nAvg AUC for {tfdir.name}: {mean_auc:.3f}")

        # Save per-TF AUC CSV
        metrics_dir = RESULTS_ROOT / tfdir.name / "pwm" / "metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)
        auc_csv_path = metrics_dir / f"{tfdir.name}_auc.csv"

        with open(auc_csv_path, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["TF", "Fold", "AUC"])
            writer.writeheader()
            writer.writerows(auc_records)
        
        print(f"AUC metrics saved to: {auc_csv_path}")

        # Add to master AUC list
        master_auc_records.extend(auc_records)
    
    # Write master AUC CSV
    master_dir = PROJECT / "master_results" / "pwm" / "metrics"
    master_dir.mkdir(parents=True, exist_ok=True)
    master_csv = master_dir / "master_auc.csv"

    with open(master_csv, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["TF", "Fold", "AUC"])
        writer.writeheader()
        writer.writerows(master_auc_records)
    print(f"Master AUC metrics saved to: {master_csv}")

print("\n====== PWM AUC calculation completed! ======")

    
    





