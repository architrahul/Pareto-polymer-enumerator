"""
leakage_analysis.py — Reproduce paper Section 7.3 leakage analysis.

Runs the covering-design Hilbert basis enumeration at user-chosen (k, t,
monomer mode) for the full linear cascade and the same cascade with the
first input monomer (x1) removed. The two Hilbert bases are then fed into
COFFEE for the equilibrium-concentration comparison shown in Figures 5–6
of the paper.

Usage:
    python leakage_analysis.py                # default: --n 8 --t 5 --k 25
    python leakage_analysis.py --n 7          # paper's |M|=56 case 
    python leakage_analysis.py --n 7 --t 3    # for Figure 5 (t=3 reduced set)
    python leakage_analysis.py --n 7 --only incomplete   # skip the "full" variant
"""

import argparse
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # parent src/

from hilbert_pipeline import (
    cleanup_normaliz_files,
    full_run_k,
    get_all_unique_domains,
    load_covering_blocks,
    load_monomers,
    save_polymer_vectors,
    start_input_listener,
)
from paths import EXAMPLE_TBNS_DIR, LEAKAGE_HB_DIR, LOGS_DIR


def leakage_analysis(cascade_n: int = 8, t: int = 5, k: int = 25,
                     only: str = "both"):
    """Run covering-design enumeration on cascade-n {full, first-input-removed}.

    only ∈ {"both", "full", "incomplete"} restricts which variant to run.
    """
    SAVE_DIR = LEAKAGE_HB_DIR
    all_files = {
        f"n{cascade_n}_full":       os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{cascade_n}.txt"),
        f"n{cascade_n}_incomplete": os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{cascade_n}_incomplete.txt"),
    }
    if only == "full":
        MONOMER_FILES = {k_: v for k_, v in all_files.items() if k_.endswith("_full")}
    elif only == "incomplete":
        MONOMER_FILES = {k_: v for k_, v in all_files.items() if k_.endswith("_incomplete")}
    else:
        MONOMER_FILES = all_files
    T    = t
    K    = k
    MODE = "monomer"

    os.makedirs(SAVE_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)
    cleanup_normaliz_files()
    start_input_listener()

    for label, monomer_file in MONOMER_FILES.items():
        print(f"\n{'='*70}")
        print(f"  System  : {label.upper()}  ({monomer_file})")
        print(f"  Mode    : {MODE}  |  t={T}  |  k={K}")
        print(f"{'='*70}")

        all_monomers = load_monomers(monomer_file)
        n_monomers   = len(all_monomers)
        all_domains  = get_all_unique_domains(all_monomers)
        n_domains    = len(all_domains)
        n            = n_monomers

        print(f"  Loaded {n_monomers} monomers, {n_domains} unique domain types")

        if K > n:
            print(f"  WARNING: k={K} > n={n} for '{label}' system — skipping.")
            continue

        log_file = os.path.join(
            SAVE_DIR,
            f"log_leakage_{MODE}_{label}_n{n_monomers}_t{T}_k{K}.txt"
        )

        with open(log_file, "a") as log:
            log.write(
                f"Leakage Analysis — Pareto-Optimal Polymer Enumeration\n"
                f"Started    : {datetime.now()}\n"
                f"System     : {label}  ({monomer_file})\n"
                f"Mode       : {MODE}  |  t={T}  |  k={K}\n"
                f"n_monomers={n_monomers}  n_domains={n_domains}\n"
                + "=" * 70 + "\n"
            )
            log.flush()

            try:
                blocks = load_covering_blocks(n, K, T, fallback_dp=False)
            except RuntimeError as e:
                msg = f"  [{label}] Could not load covering blocks: {e}\n"
                print(msg); log.write(msg); log.flush()
                continue

            print(f"  [{label}] Loaded {len(blocks)} covering blocks. Running full enumeration ...")

            try:
                result, _ = full_run_k(
                    K, blocks, n, all_monomers, MODE, all_domains, n_monomers, log
                )
            except KeyboardInterrupt:
                print("\n  Interrupted.")
                log.write("\nInterrupted by user.\n"); log.flush()
                cleanup_normaliz_files()
                return

            if result is None:
                msg = f"  [{label}] full_run_k returned None (possibly skipped/timed out).\n"
                print(msg); log.write(msg); log.flush()
                continue

            summary = (
                f"\n  [{label}] k={K}: "
                f"wall={result['total_wall_time']:.3f}s  "
                f"normaliz={result['total_normaliz_time']:.3f}s  "
                f"unique_vectors={result['unique_vectors']}\n"
            )
            print(summary); log.write(summary); log.flush()

            out_path = os.path.join(
                SAVE_DIR,
                f"hilbert_k{K}_t{T}_{MODE}_{label}.txt"
            )
            save_polymer_vectors(
                result["vectors"], out_path,
                n_monomers=n_monomers,
                comment=f"leakage analysis, covering mode, k={K}, t={T}, system={label}"
            )
            print(f"  Saved vectors → {out_path}")

    cleanup_normaliz_files()
    print("\n  Leakage analysis complete.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", "--leakage-n", type=int, default=8, dest="cascade_n",
                        help="Cascade depth. Default 8. Use 7 to reproduce the paper's |M|=56 figure.")
    parser.add_argument("--t", type=int, default=5,
                        help="Support bound t. Default 5. Use 3 for the paper's Figure 5.")
    parser.add_argument("--k", type=int, default=25,
                        help="Block size k. Default 25.")
    parser.add_argument("--only", choices=["both", "full", "incomplete"], default="both",
                        help="Which cascade variant to run. Default 'both'.")
    args = parser.parse_args()
    leakage_analysis(cascade_n=args.cascade_n, t=args.t, k=args.k, only=args.only)


if __name__ == "__main__":
    main()
