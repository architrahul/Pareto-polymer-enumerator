"""
leakage_compute_all.py — One-shot driver that computes the three Hilbert
basis files Figures 5/6 (Section 7.3) compare against each other, on the
cascade-n incomplete system (|M|=8n).

Default target: cascade n=7 incomplete (|M|=56), the paper's intended case.

The three polymer sets:
  1. P̂_{k=25, t=3}   →  Figure 5 reduced candidate
  2. P̂_{k=25, t=5}   →  Figure 6 reduced candidate
  3. Full P*              →  baseline for both figures

Each run is sequential (no concurrent Normaliz). Output files:
  results/leakage/hilbert_basis/hilbert_k25_t3_monomer_n{N}_incomplete.txt
  results/leakage/hilbert_basis/hilbert_k25_t5_monomer_n{N}_incomplete.txt
  results/leakage/hilbert_basis/hilbert_full_p_star_n{N}_incomplete.txt
"""

import argparse
import os
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hilbert_pipeline import (
    cleanup_normaliz_files,
    get_all_unique_domains,
    load_monomers,
    run_normaliz_on_subset,
    save_polymer_vectors,
)
from paths import EXAMPLE_TBNS_DIR, LEAKAGE_HB_DIR, LOGS_DIR

from leakage_analysis import leakage_analysis


def compute_full_p_star(cascade_n: int, label: str, monomer_file: str):
    """One Normaliz call on the entire monomer set = full Hilbert basis P*."""
    print(f"\n{'='*70}")
    print(f"  FULL P* (Full HB) on {label}  ({monomer_file})")
    print(f"{'='*70}")

    all_monomers = load_monomers(monomer_file)
    n            = len(all_monomers)
    all_domains  = get_all_unique_domains(all_monomers)
    print(f"  |M|={n}, |Σ|={len(all_domains)}")

    cleanup_normaliz_files()
    t0 = time.time()
    elapsed, raw_vectors = run_normaliz_on_subset(all_monomers)
    wall = time.time() - t0
    print(f"  Normaliz finished in {elapsed:.1f}s (wall {wall:.1f}s) — "
          f"{len(raw_vectors)} vectors")

    if not raw_vectors:
        print(f"  WARNING: no vectors returned (timeout or error). Skipping save.")
        return

    # Project augmented-system vectors (|M| + 2|Σ| coords) back to just the
    # M-coordinates (the π map of Theorem 9 / Lemma 12). Drop the zero vector.
    vector_tuples = {tuple(v[:n]) for v in raw_vectors}
    vector_tuples.discard(tuple([0] * n))

    out_path = os.path.join(LEAKAGE_HB_DIR, f"hilbert_full_p_star_{label}.txt")
    save_polymer_vectors(
        vector_tuples, out_path, n_monomers=n,
        comment=f"Full P* (Full HB), system={label}",
    )
    print(f"  Saved {len(vector_tuples)} vectors → {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7 (paper's |M|=56 case).")
    args = parser.parse_args()

    os.makedirs(LEAKAGE_HB_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)

    n         = args.cascade_n
    label     = f"n{n}_incomplete"
    monomers  = os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{n}_incomplete.txt")
    assert os.path.exists(monomers), f"missing monomer file: {monomers}"

    print(f"\n{'#'*70}\n# leakage_compute_all  cascade n={n} incomplete  started "
          f"{datetime.now().isoformat(timespec='seconds')}\n{'#'*70}")

    # --- Step 1: P̂_{k=25, t=3} (Figure 5) ---
    p3_file = os.path.join(LEAKAGE_HB_DIR,
                           f"hilbert_k25_t3_monomer_{label}.txt")
    if os.path.exists(p3_file):
        print(f"\n[STEP 1] {p3_file} already exists — SKIP")
    else:
        print(f"\n[STEP 1] covering @ k=25 t=3 for {label}")
        leakage_analysis(cascade_n=n, t=3, k=25, only="incomplete")

    # --- Step 2: P̂_{k=25, t=5} (Figure 6) ---
    p5_file = os.path.join(LEAKAGE_HB_DIR,
                           f"hilbert_k25_t5_monomer_{label}.txt")
    if os.path.exists(p5_file):
        print(f"\n[STEP 2] {p5_file} already exists — SKIP")
    else:
        print(f"\n[STEP 2] covering @ k=25 t=5 for {label}")
        leakage_analysis(cascade_n=n, t=5, k=25, only="incomplete")

    # --- Step 3: Full P* baseline ---
    full_file = os.path.join(LEAKAGE_HB_DIR, f"hilbert_full_p_star_{label}.txt")
    if os.path.exists(full_file):
        print(f"\n[STEP 3] {full_file} already exists — SKIP")
    else:
        print(f"\n[STEP 3] Full P* on {label}")
        compute_full_p_star(n, label, monomers)

    print(f"\n{'#'*70}\n# leakage_compute_all DONE "
          f"{datetime.now().isoformat(timespec='seconds')}\n{'#'*70}")


if __name__ == "__main__":
    main()
