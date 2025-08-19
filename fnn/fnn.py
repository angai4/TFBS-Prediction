#!/usr/bin/env python3
"""
Optuna hyperparameter tuning and 10-fold cross-validation for an FNN classifier
over multiple transcription factors (TFs).

Workflow per TF
---------------
1) Use fold_01 only to run an Optuna study and find good hyperparameters.
2) With the best params, train/evaluate a small feed-forward neural net on all
   10 folds (train on train_pos_TXX.fa vs train_neg_TXX.fa, test on test_pos_FXX.fa
   vs test_bg_BXX.fa).
3) Save:
   - ROC curves      → results/<TF>/fnn/figures/fnn_roc_fold_XX.png
   - Per-TF AUC CSV  → results/<TF>/fnn/metrics/<TF>_auc.csv
   - Detailed metrics→ results/<TF>/fnn/metrics/<TF>_detailed_metrics.txt
   - Best params     → results/<TF>/fnn/model/<TF>_best_model.txt

Also supports resuming, because the Optuna study is stored in a per-TF SQLite DB.

Assumptions
-----------
- Data rooted at ~/msc_project/data/<TF>/fnn_crossval/fold_XX with FASTA files:
    train_pos_TXX.fa, train_neg_TXX.fa, test_pos_FXX.fa, test_bg_BXX.fa
- FASTA sequences are 101 bp (or fixed length). One-hot encoding flattens them
  into a 4*L vector.
- CPU is fine; tensors created from NumPy go to CPU by default.

Notes
-----
- For reproducibility, you may want to set random seeds for numpy/torch, but
  this script does not enforce that (Optuna tries multiple trials anyway).
"""
import os
from pathlib import Path
import numpy as np
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from Bio import SeqIO
import optuna
from optuna.exceptions import TrialPruned
from optuna.storages import RDBStorage
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve
)
import matplotlib.pyplot as plt
import pandas as pd

# Global settings
BATCH_SIZE = 64
EPOCHS     = 25
DATA_ROOT  = Path("~/msc_project/data").expanduser()
RESULT_ROOT = Path("~/msc_project/results").expanduser()
#MASTER_METRICS = Path("~/msc_project/master_results/fnn").expanduser() / "metrics"
#MASTER_METRICS.mkdir(parents=True, exist_ok=True)

def create_full_sequence_dataset(pos_fasta_files, neg_fasta_files):
    """
    pos_fasta_files, neg_fasta_files : list of str or Path
        Paths to FASTA files containing positive and negative sequences.
    Returns
    -------
    seqs : list of str
        All the sequences from the pos and neg FASTAs.
    labels : list of int
        1 for each positive sequence, 0 for each negative.
    """
    # initialise empty list to store all DNA sequences
    seqs = []
    # initialise empty list to stor corresponding labels (0 or 1)
    labels = []

    # load positives
    # loop through each positive FASTA file
    for fp in pos_fasta_files:
        # parse each FASTA file and extract sequences
        for rec in SeqIO.parse(Path(fp).expanduser(), "fasta"):
            # convert sequence record to string and add to list
            seqs.append(str(rec.seq))
            # add label 1 for positive sequences
            labels.append(1)

    # load negatives
    # loop through each negative FASTA file
    for fn in neg_fasta_files:
        # parse each FASTA file and extract sequences
        for rec in SeqIO.parse(Path(fn).expanduser(), "fasta"):
            # convert sequence record to string and add to list
            seqs.append(str(rec.seq))
            # add label 0 for negative sequences
            labels.append(0)
    
    # return the sequences and labels
    return seqs, labels


def one_hot_encode(seq):
    # turn a dna sequence into a flat one-hot encoded array
    
    # define mapping from DNA bases to one-hot vectors
    mapping = {
        'A': [1.,0.,0.,0.],
        'C': [0.,1.,0.,0.],
        'G': [0.,0.,1.,0.],
        'T': [0.,0.,0.,1.],
        'N': [0.,0.,0.,0.],
    }
    # convert each base in sequence to one-hot vector, handling unknown bases
    arr = np.array([mapping.get(b, mapping['N']) for b in seq.upper()], dtype=np.float32)
    # flatten the 2D array into 1D for neural network input
    return arr.flatten()

class FeedForwardNN(nn.Module):
    # a fully connected neural network
    def __init__(self, input_dim: int, h1: int, h2: int, dropout: float):
        # initialise the parent class
        super().__init__()
        # define the neural network architecture
        self.net = nn.Sequential(
            # first layer: input to first hidden layer
            nn.Linear(input_dim, h1),
            # ReLU activation function
            nn.ReLU(),
            # dropout layer to prevent overfitting
            nn.Dropout(dropout),
            # second layer: first hidden layer to second hidden layer
            nn.Linear(h1, h2),
            # ReLU activation function
            nn.ReLU(),
            # dropout layer to prevent overfitting
            nn.Dropout(dropout),
            # output layer: second hidden laayer to single output
            nn.Linear(h2, 1),
        )
        # layers: input->h1->h2->single output (logit)
    def forward(self, x):
        # pass the data through the network
        return self.net(x).squeeze(1)


def define_model(trial, input_dim):
    # suggest hyperparameter for first hidden layer size (20-30, step 10)
    h1 = trial.suggest_int("n_units_l0", 20, 30, step=10)
    # suggest hyperparameter for second hidden layer size (5-10, step 5)
    h2 = trial.suggest_int("n_units_l1", 5, 10, step=5)
    # suggest dropout probability
    dropout = trial.suggest_float("dropout", 0.0, 0.5)
    # create and return the model with the suggested architecture
    return FeedForwardNN(input_dim, h1, h2, dropout)


def objective(trial, tf_name):
    # load fold_01 files for hyperparameter tuning 
    # define paths to positive and negative training/testing files
    base = DATA_ROOT / tf_name / "fnn_crossval" / "fold_01"
    pos = [base / "train_pos_T01.fa", base / "test_pos_F01.fa"]
    neg = [base / "train_neg_T01.fa", base / "test_bg_B01.fa"]
    # load sequenes and labels from fasta files
    seqs, labs = create_full_sequence_dataset(pos, neg)
    # convert sequences to one-hot encoded arrays and stack them
    X = np.stack([one_hot_encode(s) for s in seqs])
    # convert labels to numpy arrays with float32 type
    Y = np.array(labs, dtype=np.float32)
    # create PyTorch dataset from features and labels
    ds = TensorDataset(torch.from_numpy(X), torch.from_numpy(Y))
    # calculate size for training set (80% of data)
    n_train = int(0.8 * len(ds))
    # split dataset into training and testing sets
    train_ds, test_ds = torch.utils.data.random_split(ds, [n_train, len(ds)-n_train])

    # get input dimension from feature matrix
    input_dim = X.shape[1]
    # create model with trial-suggested hyperparameters
    model = define_model(trial, input_dim)
    # suggest learning rate hyperparameter (log scale between 1e-5 and 1e-1)
    lr    = trial.suggest_float("lr", 1e-5, 1e-1, log=True)
    # create Adam optimizer with suggested learning rate
    opt   = torch.optim.Adam(model.parameters(), lr=lr)
    # define binary cross-entropy loss function with logits
    crit  = nn.BCEWithLogitsLoss()

    # create data loader for training set with shuffling
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    # create data loader for test set without shuffling
    test_loader  = DataLoader(test_ds,  batch_size=BATCH_SIZE)

    # initialize variable to track best accuracy
    best_acc = 0.0
    # loop through each training epoch
    for _ in range(EPOCHS):
        # set model to training mode
        model.train()
        # loop through each batch in training data
        for xb, yb in train_loader:
            # calculate loss for current batch
            loss = crit(model(xb), yb)
            # clear gradients from previous step
            opt.zero_grad()
            # calculate gradients via backpropagation
            loss.backward()
            # update model parameters
            opt.step()
        
        # set model to evaluation mode
        model.eval()
        # initialise counter for correct predictions
        correct=0
        # disable gradient calculation for evaluation
        with torch.no_grad():
            # loop through each batch in test data
            for xb,yb in test_loader:
                # get predictions (sigmoid > 0.5 threshold)
                preds = (torch.sigmoid(model(xb))>0.5).float()
                # count correct predictions
                correct += (preds==yb).sum().item()
        # calculate accuracy for this epoch   
        acc = correct / len(test_loader.dataset)
        
        # report intermediate result to Optuna
        trial.report(acc, _)
        # check if trial should be pruned (early stopping)
        if trial.should_prune(): 
            raise TrialPruned()
        # update best accuracy if current is better
        best_acc = max(best_acc, acc)
    
    # return the best accuracy achieved
    return best_acc


def detailed_evaluate(best_params, tf_name):
    """
    Run 10-fold CV using the best_trial.params.
    Prints per-fold and average accuracy, precision, recall, F1, and ROC AUC.
    """
    # prepare output directories
    figs_dir = RESULT_ROOT / tf_name / "fnn/figures"
    met_dir  = RESULT_ROOT / tf_name / "fnn/metrics"
    mod_dir  = RESULT_ROOT / tf_name / "fnn/model"
    for d in (figs_dir, met_dir, mod_dir): 
        d.mkdir(parents=True, exist_ok=True)

    # for AUC csv
    auc_rows = []
    # for detailed metrics txt
    detailed_lines = []
    acc_list = []
    prec_list = []
    rec_list = []
    f1_list  = []
    auc_list = []   
    
    # initialize input dimension (will be set from first fold)
    input_dim=None

    # loop through each of the 10 folds
    for fold in range(1,11):
        base=DATA_ROOT/ tf_name/"fnn_crossval"/f"fold_{fold:02d}"
        # define paths to training and testing files for this fold
        tr_p, tr_n = base/f"train_pos_T{fold:02d}.fa", base/f"train_neg_T{fold:02d}.fa"
        te_p, te_n = base/f"test_pos_F{fold:02d}.fa" , base/f"test_bg_B{fold:02d}.fa"
        # load training and testing sequences and their labels
        seqs_tr, y_tr = create_full_sequence_dataset([tr_p], [tr_n])
        seqs_te, y_te = create_full_sequence_dataset([te_p], [te_n])
        # convert sequences to one-hot encoded arrays
        Xtr, Xte = np.stack([one_hot_encode(s) for s in seqs_tr]), np.stack([one_hot_encode(s) for s in seqs_te])
        # set input dimension from first fold's training data
        if input_dim is None: 
            input_dim=Xtr.shape[1]

        # create training/test dataset from features and labels 
        tr_ds = TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(np.array(y_tr, dtype = np.float32)))
        te_ds = TensorDataset(torch.from_numpy(Xte), torch.from_numpy(np.array(y_te, dtype = np.float32)))
        # create data loaders for training and testing datasets
        tr_ld, te_ld = DataLoader(tr_ds, BATCH_SIZE, shuffle=True), DataLoader(te_ds, BATCH_SIZE)
        
        # create model using best hyperparameters found by Optuna
        model=FeedForwardNN(
            input_dim,
            best_params['n_units_l0'],
            best_params['n_units_l1'],
            best_params['dropout']
            )
        # create optimiser with best learning rate
        opt = torch.optim.Adam(model.parameters(),lr=best_params['lr'])
        # define loss function
        crit = nn.BCEWithLogitsLoss()
        # set model to training mode
        model.train()
        # train for specified number of epochs 
        for _ in range(EPOCHS):
            # loop through training batches
            for xb, yb in tr_ld:
                # calculate loss
                loss = crit(model(xb), yb)
                # clear gradients
                opt.zero_grad()
                # backpropagate
                loss.backward()
                # update weights
                opt.step()
        
        # set model to evaluation mode
        model.eval()
        # initialize lists to store true labels, predictions, and probabilities
        y_true, y_pred, y_prob = [], [], []
        # disable gradient calculation for evaluation
        with torch.no_grad():
            # loop through test batches
            for xb, yb in te_ld:
                # get raw model outputs (logits)
                logits = model(xb)
                # convert logits to probabilities using sigmoid
                probs  = torch.sigmoid(logits)
                # convert probabilities to binary predictions (>0.5 threshold)
                preds  = (probs > 0.5).float()
                # store true labels
                y_true.extend(yb.tolist())
                # store binary predictions
                y_pred.extend(preds.tolist())
                # store probability scores
                y_prob.extend(probs.tolist())

        # compute metrics
        acc     = accuracy_score(y_true, y_pred)
        prec    = precision_score(y_true, y_pred, zero_division=0)
        rec     = recall_score(y_true, y_pred, zero_division=0)
        f1      = f1_score(y_true, y_pred, zero_division=0)
        auc     = roc_auc_score(y_true, y_prob)

        # collect for CSV
        auc_rows.append({
            "TF":    tf_name,
            "Fold":  f"fold_{fold:02d}",
            "AUC":   round(auc, 3),
        })
        auc_list.append(auc)        

        # collect for detailed metrics txt
        detailed_lines.append(
            f"Fold {fold:02d} — Acc: {acc:.3f}, Prec: {prec:.3f}, Rec: {rec:.3f}, F1: {f1:.3f}"
        )
        acc_list.append(acc)
        prec_list.append(prec)
        rec_list.append(rec)
        f1_list.append(f1)       

        # ROC plot
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        plt.figure()
        plt.plot(fpr, tpr, lw=2, label=f'ROC curve (AUC = {auc:.3f})')
        plt.plot([0,1], [0,1], color='black', linestyle='--', label='Random Classifier')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')        
        plt.title(f"{tf_name} - fold_{fold:02d} FNN ROC")
        plt.legend(loc='lower right')
        plt.grid(True)        
        plt.tight_layout()
        plt.savefig(figs_dir/f"fnn_roc_fold_{fold:02d}.png")
        plt.close()

    
    # add average row to AUC CSV
    auc_rows.append({
        "TF":   tf_name,
        "Fold": "average",
        "AUC":  round(np.mean(auc_list), 3),
    })
    df_auc = pd.DataFrame(auc_rows, columns=["TF","Fold","AUC"])
    df_auc.to_csv(met_dir / f"{tf_name}_auc.csv", index=False)

    # write detailed metrics TXT (only acc, prec, rec, f1)
    with open(met_dir / f"{tf_name}_detailed_metrics.txt", "w") as f:
        for acc, prec, rec, f1, auc, fold_line in zip(acc_list, prec_list, rec_list, f1_list, auc_list, detailed_lines):
            f.write(f"{fold_line}, AUC: {auc:.3f}\n")
        f.write("\nAverage across 10 folds:\n")
        f.write(f"  Accuracy : {np.mean(acc_list):.3f} ± {np.std(acc_list):.3f}\n")
        f.write(f"  Precision: {np.mean(prec_list):.3f} ± {np.std(prec_list):.3f}\n")
        f.write(f"  Recall   : {np.mean(rec_list):.3f} ± {np.std(rec_list):.3f}\n")
        f.write(f"  F1       : {np.mean(f1_list):.3f} ± {np.std(f1_list):.3f}\n")
        f.write(f"  AUC      : {np.mean(auc_list):.3f} ± {np.std(auc_list):.3f}\n")
            
    # save best‐model params
    with open(mod_dir / f"{tf_name}_best_model.txt", "w") as f:
        f.write("Params:\n")
        for k, v in best_params.items():
            f.write(f"  {k}: {v}\n")

    return df_auc

# Main
if __name__ == "__main__":
    #master_dfs = []

    def is_tf_already_processed(tf_name):
        """Check if TF has already been processed by looking for key output files"""
        result_dir = RESULT_ROOT / tf_name / "fnn"
        
        # Check for the main output files that indicate successful completion
        auc_file = result_dir / "metrics" / f"{tf_name}_auc.csv"
        detailed_file = result_dir / "metrics" / f"{tf_name}_detailed_metrics.txt"
        model_file = result_dir / "model" / f"{tf_name}_best_model.txt"
        
        # Also check for at least some ROC figures (indicating fold processing completed)
        figures_dir = result_dir / "figures"
        roc_files = list(figures_dir.glob("fnn_roc_fold_*.png")) if figures_dir.exists() else []
        
        # Consider processed if all key files exist and we have ROC plots
        files_exist = auc_file.exists() and detailed_file.exists() and model_file.exists()
        has_roc_plots = len(roc_files) > 0
        
        return files_exist and has_roc_plots
    
    # Count total and already processed TFs
    all_tf_dirs = [d for d in sorted(DATA_ROOT.iterdir()) if d.is_dir()]
    already_processed = [tf_dir.name for tf_dir in all_tf_dirs if is_tf_already_processed(tf_dir.name)]

    for tf_dir in all_tf_dirs:
        tf_name = tf_dir.name
        
        # Skip if already processed
        if is_tf_already_processed(tf_name):
            continue

    #for tf_dir in sorted(DATA_ROOT.iterdir()):
        #if not tf_dir.is_dir(): 
            #continue
        #tf_name = tf_dir.name
        #print(f"====== Processing {tf_name} ======")
        # study per TF
        sqlite= Path("~/msc_project/databases").expanduser()/f"optuna_fnn_{tf_name}.db"
        sqlite.parent.mkdir(exist_ok=True)
        storage = RDBStorage(
            url=f"sqlite:///{sqlite}",
            engine_kwargs={"connect_args":{"timeout":60,"check_same_thread":False}}
            )
        study = optuna.create_study(
            study_name=f"fnn_{tf_name}",
            storage=storage,
            direction="maximize",
            load_if_exists=True
            )
        
        # hyperparameter search
        study.optimize(lambda t: objective(t, tf_name),n_trials=25)
        #print("Best for", tf_name, study.best_trial.params)
        df = detailed_evaluate(study.best_trial.params, tf_name)
        #master_dfs.append(df)
    # combine all
    #master_df = pd.concat(master_dfs, ignore_index=True)
    #master_df.to_csv(MASTER_METRICS/f"master_auc.csv",index=False)
    #print("All TFs processed. Master CSV written.")
