#!/usr/bin/env python3
"""
Create 10-fold cross-validation splits for each TF in a project.

For every transcription factor (TF) directory under the provided data root:
  - Read crossval/remaining_fasta.fa (the sequences not used for top-N peaks).
  - Randomly shuffle and split into 10 folds (reproducible with random_state).
  - Write per-fold FASTAs:
        crossval/fold_XX/train_TXX.fa
        crossval/fold_XX/test_FXX.fa
    where XX ∈ {01..10} and train/test are complementary.

Inputs (command line):
  DATA_ROOT : path to the project's data directory
              e.g., ~/msc_project/data

Folder layout expected per TF:
  <DATA_ROOT>/<TF>/crossval/remaining_fasta.fa

Outputs per TF:
  <DATA_ROOT>/<TF>/crossval/fold_01/train_T01.fa
  <DATA_ROOT>/<TF>/crossval/fold_01/test_F01.fa
  ...
  <DATA_ROOT>/<TF>/crossval/fold_10/train_T10.fa
  <DATA_ROOT>/<TF>/crossval/fold_10/test_F10.fa

Dependencies:
  - Biopython (Bio.SeqIO)
  - scikit-learn (sklearn.model_selection.KFold)
"""

import os, sys                               # handles file and directory operations, and command-line arguments
from sklearn.model_selection import KFold    # provides KFold cross-validation functionality
from Bio import SeqIO                        # provides tools for reading and writing sequence data in FASTA format

if len(sys.argv) != 2:                                                      # check if the correct number of arguments is provided
    print(f"Usage: {sys.argv[0]} /path/to/project/data", file=sys.stderr)   # print usage message if not
    sys.exit(1)                                                             # exit if the argument count is incorrect

DATA_ROOT = sys.argv[1]                                                     # the root directory containing the transcription factor data
FOLDS     = 10                                                              # number of folds for cross-validation

for tf in sorted(os.listdir(DATA_ROOT)):                                    # iterate over each transcription factor directory in the data root, e.g. ~/msc_project/data
    tfdir = os.path.join(DATA_ROOT, tf)                                     # path to the transcription factor directory
    remfa = os.path.join(tfdir, "crossval", "remaining_fasta.fa")           # path to the remaining_fasta.fa file for the transcription factor
    cvdir = os.path.join(tfdir, "crossval")                                 # path to the cross-validation directory for the transcription factor

    if not os.path.isfile(remfa):                                           # check if the remaining_fasta.fa file exists
        print(f"Skipping {tf!r}: no remaining_fasta.fa")                    # if not, print skip this transcription factor
        continue                                                            # continue to the next transcription factor

    print(f"\n{tf}: loading {remfa}")                                        
    records = list(SeqIO.parse(remfa, "fasta"))                             # parse the remaining_fasta.fa file and load all sequences into a list
    n = len(records)                                                        # count the number of sequences loaded
    print(f"    -> {n} seqs")                                               # print the number of sequences loaded

    # Create fold directories directly under crossval/
    os.makedirs(cvdir, exist_ok=True)

    kf = KFold(n_splits=FOLDS, shuffle=True, random_state=42)                   # create a KFold object for 10-fold cross-validation with shuffling and a fixed random state
    # kf.split(records) = tells KFold to split 'records' (list of sequences from remaining_fasta.fa) into 10 folds
    # train_idx = indices of the training samples (i.e. which items from 'records' list to use for training)
    # test_idx = indices of the test samples
    for fold, (train_idx, test_idx) in enumerate(kf.split(records), start=1):   # loop through each fold 01-10
        fold_dir = os.path.join(cvdir, f"fold_{fold:02d}")                      # create fold directories (fold_01,...,fold_10) inside crossval/
        os.makedirs(fold_dir, exist_ok=True)         

        train_fa = os.path.join(fold_dir, f"train_T{fold:02d}.fa")              # train_fa is the path to the train_TXX.fa file                
        test_fa  = os.path.join(fold_dir, f"test_F{fold:02d}.fa")               # test_fa is the path to the test_FXX.fa file

        SeqIO.write((records[i] for i in train_idx), train_fa, "fasta")         # SeqIO.write: writes sequence records to a file, records[i] for i in train_idx: iterates over the indices in the list and selects corresponding elements from 'records' list
        SeqIO.write((records[i] for i in test_idx),  test_fa,  "fasta")

        print(f"    -> fold {fold:02d}: train={len(train_idx)} test={len(test_idx)}")

    print(f"{tf!r} done. Folds under {cvdir}")


