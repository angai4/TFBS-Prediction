# TFBS-Prediction — PWM, TFFM & FNN Evaluation

This repository contains the pipelines used for my MSc research project titled "Prediction of transcription factor binding sites: beyond position weight matrices. Using hidden Markov models and neural networks to identify motifs and predict bindings" associated with University College London. In summary the pipelines, prepare ENCODE TF ChIP-seq peak data and evaluate TFBS prediction models:

- **PWM** (Position Weight Matrices): MEME discovery, MAST scanning, GC-matched backgrounds, ROC/AUC.
- **TFFM** (Transcription Factor Flexible Models): first-order & detailed models initialized from MEME motifs.
- **FNN** (Feed-Forward Neural Network): one-hot windows with GC-matched negatives, Optuna tuning + 10-fold CV.

> **Scope note:** PWM, TFFM, and FNN have **separate scripts** and **separate Conda environments**.  
> TFFM specifically requires **Python 2.7**, **Biopython ≥ 1.61**, **GHMM**, and access to the local **TFFM** package via `PYTHONPATH`.

---

## Repository Layout

```
.
├─ pwm/                      # PWM-only scripts
├─ tffm/                     # TFFM (first-order & detailed) scripts
├─ fnn/                      # FNN preparation, tuning, and evaluation scripts
├─ envs/
│  ├─ environment.pwm.yml    # Conda env (Python 3) for PWM
│  ├─ environment.tffm.yml   # Conda env (Python 2.7) for TFFM
│  └─ environment.fnn.yml    # Conda env (Python 3) for FNN
├─ README.md
└─ .gitignore
```

---

## Data & Reference Assumptions

Scripts assume the following under your home directory:

```
~/msc_project/
├─ data/<TF>/
│  ├─ raw/                   # downloaded peak files (ENCODE)
│  ├─ sorted/                # concatenated + sorted peaks
│  ├─ top600/                # ±50bp windows around top 600 summits (+ FASTA)
│  ├─ crossval/              # 10-fold CV over remaining peaks
│  └─ fnn_crossval/          # FNN-friendly copies (train/test pos/neg)
├─ results/<TF>/...          # per-method outputs (motif/tffm_*/fnn)
├─ master_results/...        # cross-TF master metrics (CSV)
├─ background/...            # GC-unique windows, GC% catalogs, etc.
└─ genomes/hg19.fa           # reference FASTA (see below)
```

**Reference genome path** expected by many scripts:

```bash
mkdir -p ~/genomes
# place your hg19 FASTA here:
#   ~/genomes/hg19.fa
samtools faidx ~/genomes/hg19.fa
```

---

## Environments (Conda)

### A) PWM environment (Python 3)

```bash
conda env create -f envs/environment.pwm.yml
conda activate tfbs-pred-pwm
```

Verify:

```bash
python -c "import sys, Bio, sklearn; print(sys.version); print('biopython:', Bio.__version__, 'sklearn:', sklearn.__version__)"
bedtools --version
samtools --version | head -n1
meme -version
mast -version
bigWigToBedGraph 2>&1 | head -n1
jq --version
```

`envs/environment.pwm.yml` pins common CLI tools (bedtools, samtools, MEME+MAST, ucsc-bigwigtobedgraph, jq, etc.).

### B) TFFM environment (Python 2.7)

Requirements: Python 2.7, Biopython ≥ 1.61 (Py2-compatible), the GHMM library (`import ghmm`), and the TFFM package on your `PYTHONPATH`.

```bash
conda env create -f envs/environment.tffm.yml
conda activate tfbs-pred-tffm

# Point Python to your local TFFM package (adjust path if different)
export PYTHONPATH="/home/aaron/msc_project/TFFM:${PYTHONPATH}"
```

Quick checks:

```bash
python2 -c "import sys; print(sys.version)"
python2 -c "import Bio; print('biopython:', Bio.__version__)"
python2 -c "import ghmm; print('ghmm ok')"
python2 -c "import tffm_module, constants; print('TFFM ok')"
```

If `pip install ghmm` fails on your platform, install GHMM via your OS or build from source; then re-try `import ghmm`.

### C) FNN environment (Python 3)

```bash
conda env create -f envs/environment.fnn.yml
conda activate tfbs-pred-fnn
```

Includes: numpy, scipy, pandas, matplotlib, scikit-learn, biopython, pytorch (CPU), optuna.

---

## PWM Evaluation Workflow 

1. **Download & prepare peaks**

```bash
bash pwm/encode_extraction.sh
bash pwm/cat_sort.sh
bash pwm/filter_lt_1800.sh         # drops TFs with < 1800 peaks
```

2. **Create FASTAs**

```bash
bash pwm/bed_to_fasta.sh       # ±50 bp around top 600 summits → MEME input
bash pwm/remaining_bed_to_fasta.sh    # remaining peaks → crossval/remaining_fasta.fa
```

3. **10-fold cross-validation splits**

```bash
python3 pwm/make_cv_kfold.py ~/msc_project/data
```

4. **Motif discovery & scanning**

```bash
bash pwm/meme_v2.sh         # MEME on top600 FASTA per TF
bash pwm/run_mast.sh          # MAST on train sets per fold
```

5. **Background windows & PWM evaluation**

```bash
python3 pwm/calc_gc_crg36.py
python3 pwm/make_test_bed.py
python3 pwm/calc_test_gc.py
python3 pwm/make_bg_fa_v2.py
python3 pwm/pwm_auc_calc.py              # per-fold ROC/AUC, per-TF & master CSVs
```

**Outputs:**
- Per TF: `results/<TF>/motif/…`, `results/<TF>/pwm/{figures,metrics}`
- Master: `master_results/pwm/metrics/master_auc.csv`

---

## TFFM Evaluation Workflow

Activate TFFM env and set PYTHONPATH (see above), then:

### First-order model

```bash
python2 tffm/tffm_1st_order_v2.4.py \
  [--max-tfs N] \
  [--tfs TF1 TF2 ...]
```

### Detailed model

```bash
python2 tffm/TFFM_detailed.py \
  [--max-tfs N] \
  [--tfs TF1 TF2 ...]
```

**Inputs:**
```
data/<TF>/crossval/fold_XX/{train_TXX.fa, test_FXX.fa, bg_BXX.fa}
results/<TF>/motif/meme.meme
```

**Outputs:**
- First-order: `results/<TF>/tffm_1st_order/{model,figures,metrics}`
- Detailed: `results/<TF>/tffm_detailed/{model,figures,metrics}`
- Master CSVs at:
  - `master_results/tffm_1st_order/metrics/master_auc.csv`
  - `master_results/tffm_detailed/metrics/master_auc.csv`

---

## FNN Evaluation Workflow

1. **Prepare FNN CV copies** (renames CV files into FNN-friendly names)

```bash
bash fnn/create_fnn_crossval.sh
```

2. **Create BED from training positives** (parse FASTA headers)

```bash
python3 fnn/make_train_bed.py ~/msc_project/data
```

3. **Compute GC% for training positives**

```bash
python3 fnn/calc_train_gc.py
# writes: data/<TF>/fnn_crossval/fold_XX/train_pos_TXX.gc.txt
```

4. **GC-matched training negatives**

```bash
python3 fnn/make_train_neg_fa.py
# writes: data/<TF>/fnn_crossval/fold_XX/train_neg_TXX.fa
```

5. **Optuna tuning (fold_01) + Full 10-fold CV**

```bash
python3 fnn/fnn.py
# per TF:
#   results/<TF>/fnn/figures/fnn_roc_fold_*.png
#   results/<TF>/fnn/metrics/<TF>_auc.csv
#   results/<TF>/fnn/metrics/<TF>_detailed_metrics.txt
#   results/<TF>/fnn/model/<TF>_best_model.txt
```

6. **Aggregate across TFs**

```bash
python3 fnn/write_fnn_master_csv.py
# writes: master_results/fnn/metrics/master_auc.csv
```

---

## Troubleshooting

- **"command not found" for bedtools/samtools/meme/mast**  
  Ensure you're in the PWM env: `conda activate tfbs-pred-pwm`. Then:
  ```bash
  conda install -c conda-forge -c bioconda bedtools samtools meme ucsc-bigwigtobedgraph jq
  ```

- **TFFM import fails (tffm_module/constants not found)**  
  Set your `PYTHONPATH` to the local TFFM source:
  ```bash
  export PYTHONPATH="/path/to/TFFM:${PYTHONPATH}"
  ```

- **ModuleNotFoundError: ghmm in TFFM env**  
  Try:
  ```bash
  conda activate tfbs-pred-tffm
  pip install ghmm  # or build/install OS package; ensure `import ghmm` works
  ```

- **Matplotlib display errors**  
  TFFM & evaluation scripts set `Agg` backend; no display is required.

- **Different genome path**  
  Update script constants (`GENOME`, `CRG_ROOT`, etc.) to match your setup.
