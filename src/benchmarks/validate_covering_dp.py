"""
validate_covering_dp.py — Validate compute_covering_dp against reference
covering designs in covering_design_validation/.

For each `C{v}_{k}_{t}_actual_cover.csv` file:
  1. Load the reference blocks; verify the file is itself a valid (v,k,t)-covering.
  2. Run compute_covering_dp(v, k, t) and verify its output is a valid
     (v,k,t)-covering.
  3. Print sizes side-by-side.

Note: these (v, k, t) parameters have k > 25 which is outside the LJCR
range, so the implementation falls through to the GPK Section 5
construction — exactly the code we want to validate.
"""

import csv
import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hilbert_pipeline import compute_covering_dp, _COVERING_DP_CACHE
from paths import REPO_ROOT

VAL_DIR = os.path.join(REPO_ROOT, "covering_design_validation")


def load_reference(csv_path):
    """Return (v, k, t, blocks) where blocks is a list of int lists."""
    blocks = []
    v = k = t = None
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            v = int(row["v"]); k = int(row["k"]); t = int(row["t"])
            blocks.append([int(x) for x in row["block"].split()])
    return v, k, t, blocks


def is_valid_covering(blocks, v, k, t):
    """Verify every t-subset of [1..v] is contained in some block."""
    block_sets = [frozenset(b) for b in blocks]
    for T in itertools.combinations(range(1, v + 1), t):
        Tset = set(T)
        if not any(Tset <= B for B in block_sets):
            return False, T
    return True, None


def main():
    files = sorted(f for f in os.listdir(VAL_DIR) if f.endswith(".csv"))
    print(f"Validating {len(files)} reference files in {VAL_DIR}\n")
    print(f"{'file':<35} {'(v,k,t)':<12} {'ref_size':>9} {'ref_valid':>10} "
          f"{'dp_size':>9} {'dp_valid':>10} {'dp_secs':>9}")
    print("-" * 100)

    for fname in files:
        path = os.path.join(VAL_DIR, fname)
        v, k, t, ref_blocks = load_reference(path)
        ref_valid, _ = is_valid_covering(ref_blocks, v, k, t)

        _COVERING_DP_CACHE.clear()
        t0 = time.time()
        try:
            dp_blocks = compute_covering_dp(v, k, t)
            elapsed = time.time() - t0
            dp_valid, _ = is_valid_covering(dp_blocks, v, k, t)
            dp_size = len(dp_blocks)
        except Exception as e:
            elapsed = time.time() - t0
            print(f"{fname:<35} ({v},{k},{t})  ERROR: {e}")
            continue

        print(f"{fname:<35} ({v},{k},{t}) "
              f"{len(ref_blocks):>9d} {str(ref_valid):>10} "
              f"{dp_size:>9d} {str(dp_valid):>10} {elapsed:>8.1f}s")


if __name__ == "__main__":
    main()
