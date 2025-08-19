# TFBS-Prediction — PWM Evaluation

> **Scope:** This repository currently contains pipelines and utilities for **PWM-based** motif discovery, scanning, and evaluation.  
> **Planned:** Separate modules/environments for **HMM/TFFM** and **FNN** will be added later.

---

## Setup Instructions (Conda)

### Create and activate the environment
Use the provided environment file:

```bash
conda env create -f envs/environment.pwm.yml
conda activate tfbs-pred-pwm
python -c "import sys, Bio, sklearn; print(sys.version); print('biopython:', Bio.__version__, 'sklearn:', sklearn.__version__)"
bedtools --version
samtools --version | head -n1
meme -version
mast -version
bigWigToBedGraph 2>&1 | head -n1
jq --version
mkdir -p ~/genomes
samtools faidx ~/genomes/hg19.fa
PWM Workflow (high level)

Download & prepare peaks
download_encode_peaks.sh → concat_sort_peaks.sh → prune_low_peaks.sh

Create FASTAs
make_top600_fasta.sh (for MEME)
make_remaining_fasta.sh (rest for CV)

Cross-validation splits
make_folds.py

Motif discovery & scanning
run_meme_top600.sh → run_mast_train.sh

Backgrounds & evaluation
build_background_windows.py → gc_from_folds.py → gc_matched_background.py → pwm_auc.py

Roadmap

envs/environment.hmm.yml + scripts for HMM/TFFM

envs/environment.fnn.yml + scripts for FNN
