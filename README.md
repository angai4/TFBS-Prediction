# TFBS-Prediction — PWM & TFFM Evaluation

This repository contains pipelines to prepare ENCODE TF peak data and evaluate motif models.

-  **Current modules**
  - **PWM** (Position Weight Matrices): discovery (MEME), scanning (MAST), background matching, ROC/AUC evaluation.
  - **TFFM** (Transcription Factor Flexible Models): first-order and detailed models initialized from MEME motifs.

-  **Planned**: separate modules/environments for HMM/TFFM extensions and FNN baselines.

---

## Repository Layout

```
.
├─ pwm/                      # PWM-only scripts
├─ tffm/                     # TFFM (first-order & detailed) scripts
├─ envs/
│  ├─ environment.pwm.yml    # Conda env for PWM (Python 3)
│  └─ environment.tffm.yml   # Conda env for TFFM (Python 2.7)
├─ README.md
└─ .gitignore
```

---

## Setup Instructions (Conda)

### A) PWM environment (Python 3)

1) Create the environment:
```bash
conda env create -f envs/environment.pwm.yml
conda activate tfbs-pred-pwm
```

2) Verify tools:
```bash
python -c "import sys, Bio, sklearn; print(sys.version); print('biopython:', Bio.__version__, 'sklearn:', sklearn.__version__)"
bedtools --version
samtools --version | head -n1
meme -version
mast -version
bigWigToBedGraph 2>&1 | head -n1
jq --version
```

3) Reference genome assumed by scripts:
```bash
mkdir -p ~/genomes
# place hg19.fa as ~/genomes/hg19.fa
samtools faidx ~/genomes/hg19.fa
```

**envs/environment.pwm.yml**

```yaml
name: tfbs-pred-pwm
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - numpy
  - scipy
  - matplotlib
  - scikit-learn
  - biopython
  - pip
  - bedtools>=2.31
  - samtools
  - meme=5.*
  - ucsc-bigwigtobedgraph
  - jq
  - wget
  - curl
  - gawk
  - coreutils
```

### B) TFFM environment (Python 2.7)

**Requirements:**
- Python 2.7
- Biopython ≥ 1.61
- GHMM library (Python module `ghmm` importable)
- TFFM package available on your `PYTHONPATH` (e.g., the folder that contains `tffm_module.py`)

1) Create & activate the TFFM env (legacy packages):
```bash
conda env create -f envs/environment.tffm.yml
conda activate tfbs-pred-tffm
```

2) Point Python to your TFFM package:
```bash
# Example: if your TFFM source lives at /home/aaron/msc_project/TFFM
export PYTHONPATH="/home/aaron/msc_project/TFFM:${PYTHONPATH}"
```
*The scripts also prepend that path with `sys.path.insert(...)`; update there if your path differs.*

3) Verify imports:
```bash
python2 -c "import sys; print(sys.version)"
python2 -c "import Bio;  import numpy, sklearn, matplotlib; print('Bio:', Bio.__version__)"
python2 -c "import ghmm; print('ghmm ok')"
python2 -c "import tffm_module, constants; print('TFFM ok')"
```

**envs/environment.tffm.yml** *(versions chosen to match the last Python-2-compatible releases)*:

```yaml
name: tfbs-pred-tffm
channels:
  - defaults
  - conda-forge
dependencies:
  - python=2.7
  - numpy=1.16.*
  - scipy=1.2.*
  - matplotlib=2.2.*
  - scikit-learn=0.20.4
  - biopython>=1.61,<1.78
  - pip
  # GHMM is often not available via conda; try pip or your OS package:
  #   pip install ghmm
  # or build from source following GHMM instructions.
```

 **Notes for TFFM/GHMM on Python 2.7**
- Some mirrors no longer host Py2 packages. If conda fails, try pip within the env.
- GHMM may require system dependencies (e.g., `libghmm`) and compilation. Ensure `import ghmm` works before running TFFM scripts.

---

## Data Assumptions

Scripts expect the following structure under `~/msc_project`:

```bash
data/<TF>/crossval/fold_XX/
  ├─ train_TXX.fa
  ├─ test_FXX.fa
  └─ bg_BXX.fa

results/<TF>/motif/meme.meme

genomes/hg19.fa  (indexed: hg19.fa.fai)
```

---

## Usage

### PWM workflow (high level)

1. **Download & prepare peaks**
   `pwm/download_encode_peaks.sh` → `pwm/concat_sort_peaks.sh` → `pwm/prune_low_peaks.sh`

2. **Create FASTAs**
   - `pwm/make_top600_fasta.sh` (MEME input)
   - `pwm/make_remaining_fasta.sh` (rest for CV)

3. **Cross-validation splits**
   `pwm/make_folds.py`

4. **Motif discovery & scanning**
   `pwm/run_meme_top600.sh` → `pwm/run_mast_train.sh`

5. **Backgrounds & evaluation**
   `pwm/build_background_windows.py` → `pwm/gc_from_folds.py` → `pwm/gc_matched_background.py` → `pwm/pwm_auc.py`

### TFFM evaluation

**First-order TFFM:**
```bash
# activate the Python 2.7 env and set PYTHONPATH to your TFFM package
conda activate tfbs-pred-tffm
export PYTHONPATH="/home/aaron/msc_project/TFFM:${PYTHONPATH}"

python2 tffm/tffm_first_order_cv.py \
  [--max-tfs N] \
  [--tfs TF1 TF2 ...]
```

**Detailed TFFM:**
```bash
conda activate tfbs-pred-tffm
export PYTHONPATH="/home/aaron/msc_project/TFFM:${PYTHONPATH}"

python2 tffm/tffm_detailed_cv.py \
  [--max-tfs N] \
  [--tfs TF1 TF2 ...]
```

---

## Outputs

**Per TF:**
- First-order: `results/<TF>/tffm_1st_order/{model,figures,metrics}`
- Detailed: `results/<TF>/tffm_detailed/{model,figures,metrics}`

**Master CSVs:**
- `master_results/tffm_1st_order/metrics/master_auc.csv`
- `master_results/tffm_detailed/metrics/master_auc.csv`

---

## Troubleshooting

- **ModuleNotFoundError: ghmm**  
  Install GHMM inside the `tfbs-pred-tffm` env (via `pip install ghmm` or system package), then re-run `python2 -c "import ghmm"`.

- **TFFM import fails**  
  Ensure `PYTHONPATH` includes the directory with `tffm_module.py` and `constants.py`.

- **matplotlib backend errors**  
  Scripts set `Agg` backend for headless plotting; no display is required.

- **Different genome path**  
  Update `GENOME`/`CRG_ROOT` variables in scripts to match your setup.
