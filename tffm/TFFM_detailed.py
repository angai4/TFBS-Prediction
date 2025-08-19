#!/usr/bin/env python2
"""
TFFM (Detailed) cross-validation scorer

What this script does (per TF, per fold):
  1) Initialize a Detailed (higher-order) TFFM from a MEME motif.
  2) Train the TFFM on fold-specific training FASTA.
  3) Scan positives (test) and negatives (GC-matched background) and extract
     per-sequence probabilities.
  4) Compute ROC and AUC; save a ROC plot per fold.
  5) Save per-TF AUC CSVs (fold rows only) and a master CSV across all TFs.

Assumed project layout under ~/msc_project:
  data/<TF>/crossval/fold_XX/{train_TXX.fa, test_FXX.fa, bg_BXX.fa}
  results/<TF>/motif/meme.meme
  results/<TF>/tffm_detailed/{model,figures,metrics}
  master_results/tffm_detailed/metrics/master_auc.csv

Notes:
  - Uses Python 2 (for TFFM compatibility in many setups).
  - Matplotlib backend set to 'Agg' for headless plotting.
"""

import os
import sys
import csv
import gc
import logging
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add TFFM directory to Python path (adjust as needed)
sys.path.insert(0, "/home/aaron/msc_project/TFFM")
import tffm_module
from constants import TFFM_KIND

# Set up logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )


def ensure_dir(path):
    """Create directory and parents if not exists."""
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def extract_probabilities(hits):
    """Extract probability values (last column) from TFFM hits."""
    probs = []
    for hit in hits:
        if not hit:
            continue
        parts = str(hit).split('\t')
        if len(parts) >= 8:
            try:
                probs.append(float(parts[-1]))
            except ValueError:
                pass
    return probs


def process_fold(args):
    (tf, fold_name, meme_file,
     model_dir, figures_dir,
     train_seqs, test_seqs, bg_seqs,
     fold_number) = args

    # Train initial model
    tffm = tffm_module.tffm_from_meme(meme_file, TFFM_KIND.DETAILED)
    init_xml = os.path.join(model_dir, "tffm_detailed_initial_fold_{}.xml".format(fold_number))
    tffm.write(init_xml)
    tffm.train(train_seqs)
    trained_xml = os.path.join(model_dir, "tffm_detailed_fold_{}.xml".format(fold_number))
    tffm.write(trained_xml)
    logging.info("%s %s: Model trained and saved to %s", tf, fold_name, trained_xml)
    del tffm

    # Load trained model
    tffm = tffm_module.tffm_from_xml(trained_xml, TFFM_KIND.DETAILED)

    # Scan sequences
    pos_hits = list(tffm.scan_sequences(test_seqs, only_best=True))
    neg_hits = list(tffm.scan_sequences(bg_seqs, only_best=True))
    pos_probs = extract_probabilities(pos_hits)
    neg_probs = extract_probabilities(neg_hits)
    del pos_hits, neg_hits, tffm
    gc.collect()

    # Check data
    if not pos_probs or not neg_probs:
        logging.error("%s %s: No probabilities found", tf, fold_name)
        return {"TF": tf, "Fold": fold_name, "AUC": None}

    # ROC and AUC
    y_true = np.array([1] * len(pos_probs) + [0] * len(neg_probs))
    y_scores = np.array(pos_probs + neg_probs)
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    logging.info("%s %s: AUC = %.3f", tf, fold_name, roc_auc)

    # Plot ROC
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label='ROC curve (AUC = %.3f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--', label='Random classifier')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('%s - %s TFFM Detailed ROC' % (tf, fold_name))
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'tffm_detailed_roc_%s.png' % fold_name)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info("%s %s: ROC plot saved to %s", tf, fold_name, fig_path)

    return {"TF": tf, "Fold": fold_name, "AUC": round(roc_auc, 3)}


def main():
    setup_logging()

    parser = argparse.ArgumentParser(description="Run TFFM detailed cross-val and ROC/AUC.")
    parser.add_argument(
        '--max-tfs', type=int, default=None,
        help='If set, only process this many TFs (in sorted order)'
    )
    parser.add_argument(
        '--tfs', nargs='+', default=None,
        help='If set, only process these TF names'
    )
    args = parser.parse_args()

    # Directories
    home = os.path.expanduser('~')
    project = os.path.join(home, 'msc_project')
    data_root = os.path.join(project, 'data')
    results_root = os.path.join(project, 'results')
    masters_root = os.path.join(project, 'master_results', 'tffm_detailed', 'metrics')
    ensure_dir(masters_root)

    # Prepare jobs
    jobs = []
    seen = 0
    for tfdir in sorted(os.listdir(data_root)):
        # whitelist specific TFs
        if args.tfs and tfdir not in args.tfs:
            continue
        # cap number of TFs
        if args.max_tfs and seen >= args.max_tfs:
            break
        seen += 1

        tf_path = os.path.join(data_root, tfdir)
        crossval = os.path.join(tf_path, 'crossval')
        meme_file = os.path.join(results_root, tfdir, 'motif', 'meme.meme')
        if not (os.path.isdir(tf_path) and os.path.isdir(crossval) and os.path.isfile(meme_file)):
            continue
        logging.info("Scheduling TF: %s", tfdir)

        # Prepare per-TF dirs
        tffm_root = os.path.join(results_root, tfdir, 'tffm_detailed')
        model_dir = os.path.join(tffm_root, 'model')
        figures_dir = os.path.join(tffm_root, 'figures')
        metrics_dir = os.path.join(tffm_root, 'metrics')
        ensure_dir(model_dir)
        ensure_dir(figures_dir)
        ensure_dir(metrics_dir)

        for fold in sorted(os.listdir(crossval)):
            fold_dir = os.path.join(crossval, fold)
            if not os.path.isdir(fold_dir):
                continue
            fold_num = fold.split('_')[-1]
            train = os.path.join(fold_dir, 'train_T%s.fa' % fold_num)
            test  = os.path.join(fold_dir, 'test_F%s.fa' % fold_num)
            bg    = os.path.join(fold_dir, 'bg_B%s.fa'   % fold_num)
            if not (os.path.isfile(train) and os.path.isfile(test) and os.path.isfile(bg)):
                logging.warning("%s %s: Missing train/test/bg files", tfdir, fold)
                continue
            jobs.append((
                tfdir, fold, meme_file,
                model_dir, figures_dir,
                train, test, bg,
                fold_num
            ))

    # Run folds in parallel, capped to 70% of cores
    max_procs = 40
    pool_size = min(len(jobs), max_procs)
    logging.info("Creating pool with %d worker(s)", pool_size)
    pool = Pool(processes=pool_size)

    results = pool.map(process_fold, jobs)
    pool.close()
    pool.join()
    logging.info("Got %d results back from pool", len(results))

    # Group results by TF
    tf_groups = {}
    for rec in results:
        tf_groups.setdefault(rec['TF'], []).append(rec)

    master_records = []
    for tf, recs in tf_groups.items():
        valid = [r for r in recs if r['AUC'] is not None]
        if not valid:
            logging.warning("%s: No valid AUC records", tf)
            continue
        # per-TF CSV
        csv_path = os.path.join(results_root, tf, 'tffm_detailed', 'metrics', '%s_auc.csv' % tf)
        with open(csv_path, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['TF', 'Fold', 'AUC'])
            writer.writeheader()
            writer.writerows(valid)
        logging.info("%s: Wrote per-TF AUC CSV to %s", tf, csv_path)
        # average
        avg_auc = round(sum(r['AUC'] for r in valid) / len(valid), 3)
        avg_rec = {'TF': tf, 'Fold': 'Average', 'AUC': avg_auc}
        valid.append(avg_rec) # this didnt work for per-TF CSV 
        master_records.extend(valid)

    # master CSV
    master_csv = os.path.join(masters_root, 'master_auc.csv')
    with open(master_csv, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['TF', 'Fold', 'AUC'])
        writer.writeheader()
        writer.writerows(master_records)
    logging.info("Wrote master AUC CSV to %s", master_csv)

    logging.info("All processing complete.")

if __name__ == '__main__':
    main()
