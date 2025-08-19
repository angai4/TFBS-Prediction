#!/usr/bin/env python3
"""
Optimized background sequence sampler for ML training data generation.

Purpose
-------
Generate GC-content–matched *negative* sequences for TF binding prediction
experiments. For each TF and fold, the script reads the GC% distribution of
positive training windows (computed previously from BED via `bedtools nuc`)
and samples an equal number of background windows per GC% bin from a large,
precomputed pool of candidate background sequences.

High-level flow
---------------
1) Load a background catalog with entries of the form:
      chr  start  end  gc%  sequence
   and group all entries by integer GC% bin (e.g., 37 → all ~37% GC).

2) For each TF's `fnn_crossval/fold_XX/`:
   - Read `train_pos_TXX.gc.txt` to count how many positives exist in each
     GC% bin for that fold.
   - For each required bin, sample the same number of background windows from
     the background catalog. If the exact bin lacks enough entries, fallback to
     lower GC bins (target-1, target-2, …) up to `MAX_FALLBACK_DISTANCE`.

3) Write the sampled negatives to FASTA:
      fnn_crossval/fold_XX/train_neg_TXX.fa
   with headers in the form: >chr:start-end

Inputs (expected layout)
------------------------
~/msc_project/
  ├── background/windows.gc.seq.txt      # background catalog (chr, start, end, gc, seq)
  └── data/<TF>/fnn_crossval/fold_XX/
        └── train_pos_TXX.gc.txt         # per-positive GC% (chr, start, end, gc%)

Outputs
-------
~/msc_project/data/<TF>/fnn_crossval/fold_XX/train_neg_TXX.fa

Notes
-----
- GC% binning uses integer truncation of the GC column in inputs (e.g., 37.9 → 37).
- Fallback strategy only searches *downward* in GC (target-1, target-2, …).
- Set `MAX_FALLBACK_DISTANCE` to control the maximum GC% distance for fallback.
"""
import sys
import random
import logging
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional

# Configuration
DATA_ROOT = Path.home() / "msc_project" / "data"
BG_ROOT = Path.home() / "msc_project" / "background" / "windows.gc.seq.txt"
MAX_FALLBACK_DISTANCE = 50

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_background_sequences(bg_file: Path) -> Dict[int, List[Tuple[str, str, str, str]]]:
    """
    Load and group background sequences by GC-bin.
    
    Args:
        bg_file: Path to background sequences file
        
    Returns:
        Dictionary mapping GC% bins to lists of (chr, start, end, seq) tuples
    """
    logger.info(f"Loading background sequences from {bg_file}")
    bg_by_bin = defaultdict(list)
    
    try:
        with open(bg_file) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    parts = line.split(None, 4)
                    if len(parts) != 5:
                        logger.warning(f"Line {line_num}: Expected 5 fields, got {len(parts)}")
                        continue
                        
                    chr_, start, end, gc, seq = parts
                    bin_gc = int(float(gc))
                    bg_by_bin[bin_gc].append((chr_, start, end, seq))
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Line {line_num}: Parse error - {e}")
                    continue
                    
    except FileNotFoundError:
        logger.error(f"Background file not found: {bg_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading background file: {e}")
        sys.exit(1)
    
    logger.info(f"Loaded {len(bg_by_bin)} GC% bins from background data")
    
    # Log bin sizes for diagnostics
    total_seqs = sum(len(seqs) for seqs in bg_by_bin.values())
    logger.info(f"Total background sequences: {total_seqs}")
    
    return dict(bg_by_bin)  # Convert to regular dict


def find_sequences_with_fallback(
    target_gc: int, 
    needed_count: int, 
    bg_by_bin: Dict[int, List[Tuple[str, str, str, str]]], 
    max_fallback_distance: int = MAX_FALLBACK_DISTANCE
) -> Optional[List[Tuple[str, str, str, str]]]:
    """
    Find sequences for target GC%, with fallback to nearby GC% if needed.
    
    Args:
        target_gc: Target GC percentage
        needed_count: Number of sequences needed
        bg_by_bin: Dictionary of sequences by GC bin
        max_fallback_distance: Maximum GC% distance to fall back
    
    Returns:
        List of selected sequences, or None if insufficient data
    """
    for gc_offset in range(max_fallback_distance + 1):
        fallback_gc = target_gc - gc_offset
        
        if fallback_gc < 0:
            continue
            
        pool = bg_by_bin.get(fallback_gc, [])
        
        if len(pool) >= needed_count:
            if gc_offset > 0:
                logger.info(f"  Fallback: Using {fallback_gc}% GC (target was {target_gc}%)")
            return random.sample(pool, needed_count)
    
    return None


def count_gc_content(gc_file: Path) -> Counter:
    """
    Count sequences per GC-bin from a GC content file.
    
    Args:
        gc_file: Path to GC content file
        
    Returns:
        Counter object with GC% bin counts
    """
    counts = Counter()
    
    try:
        with open(gc_file) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                try:
                    parts = line.split()
                    if len(parts) < 4:
                        logger.warning(f"Line {line_num}: Expected at least 4 fields, got {len(parts)}")
                        continue
                        
                    gc = float(parts[3])
                    bin_gc = int(gc)
                    counts[bin_gc] += 1
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Line {line_num}: Parse error - {e}")
                    continue
                    
    except FileNotFoundError:
        logger.error(f"GC file not found: {gc_file}")
        return Counter()
    except Exception as e:
        logger.error(f"Error reading GC file: {e}")
        return Counter()
        
    return counts


def write_fasta_sequences(sequences: List[Tuple[str, str, str, str]], output_file: Path) -> None:
    """
    Write sequences to FASTA format file.
    
    Args:
        sequences: List of (chr, start, end, seq) tuples
        output_file: Output FASTA file path
    """
    try:
        with open(output_file, 'w') as f:
            for chr_, start, end, seq in sequences:
                f.write(f">{chr_}:{start}-{end}\n{seq}\n")
        logger.info(f"Wrote {len(sequences)} sequences to {output_file}")
    except Exception as e:
        logger.error(f"Error writing FASTA file {output_file}: {e}")
        raise


def process_tf_fold(fold_dir: Path, bg_by_bin: Dict[int, List[Tuple[str, str, str, str]]]) -> bool:
    """
    Process a single fold directory for a transcription factor.
    
    Args:
        fold_dir: Path to fold directory
        bg_by_bin: Background sequences grouped by GC%
        
    Returns:
        True if successful, False if should skip TF
    """
    fold_number = fold_dir.name.split('_')[-1]
    name = f"T{fold_number.zfill(2)}"
    gc_file = fold_dir / f"train_pos_{name}.gc.txt"
    
    logger.info(f"Processing fold: {fold_dir.name}")
    
    # Count GC content distribution
    gc_counts = count_gc_content(gc_file)
    if not gc_counts:
        logger.error(f"No valid GC data found in {gc_file}")
        return False
        
    logger.info(f"Found {len(gc_counts)} GC% bins, total sequences: {sum(gc_counts.values())}")
    
    # Sample background sequences with fallback
    selected_sequences = []
    
    for bin_gc, needed_count in gc_counts.items():
        sequences = find_sequences_with_fallback(bin_gc, needed_count, bg_by_bin)
        
        if sequences is None:
            logger.error(f"Insufficient background data for {bin_gc}% GC (need {needed_count})")
            logger.error("Skipping TF due to insufficient background data")
            return False
            
        selected_sequences.extend(sequences)
    
    # Write output FASTA
    output_file = fold_dir / f"train_neg_{name}.fa"
    write_fasta_sequences(selected_sequences, output_file)
    
    return True


def main():
    """Main processing function."""
    logger.info("Starting background sequence sampling")
    
    # Load background sequences
    bg_by_bin = load_background_sequences(BG_ROOT)
    if not bg_by_bin:
        logger.error("No background sequences loaded")
        sys.exit(1)
    
    # Process each TF
    processed_tfs = 0
    skipped_tfs = 0
    
    for tf_dir in sorted(DATA_ROOT.iterdir()):
        if not tf_dir.is_dir():
            continue
            
        fnn_crossval = tf_dir / "fnn_crossval"
        if not fnn_crossval.is_dir():
            continue
        
        logger.info(f"==== Processing TF: {tf_dir.name} ====")
        
        tf_success = True
        fold_dirs = list(fnn_crossval.glob("fold_[0-9][0-9]"))
        
        if not fold_dirs:
            logger.warning(f"No fold directories found for {tf_dir.name}")
            continue
        
        for fold_dir in sorted(fold_dirs):
            if not process_tf_fold(fold_dir, bg_by_bin):
                tf_success = False
                break
        
        if tf_success:
            processed_tfs += 1
            logger.info(f"Successfully processed TF: {tf_dir.name}")
        else:
            skipped_tfs += 1
            logger.warning(f"Skipped TF: {tf_dir.name}")
        
        logger.info("")  # Blank line for readability
    
    # Summary
    logger.info(f"Processing complete: {processed_tfs} TFs processed, {skipped_tfs} TFs skipped")


if __name__ == "__main__":
    main()