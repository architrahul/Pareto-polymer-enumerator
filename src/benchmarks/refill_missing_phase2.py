"""
refill_missing_phase2.py — Fill in the Phase-2 cells (cascade family, t)
that previously had best_k=None (no covering available) using:

    • k = 25 fixed
    • probe-only (no full enumeration) — record `best_projected` only
    • fallback_dp = True for binary d=4 (n=121, outside LJCR's n<100 range);
      LJCR-only is enough for dna m=6/7 (both n<100).

Updates results/experiments/phase2_covering.json in place: each refilled
record gets best_k=25, best_projected=<probe total>, and
full_run_skipped_reason="refilled_probe_only" so make_plots.py treats it as
a probe-projection bar (hatched / "est").
"""

import json
import os
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hilbert_pipeline import (
    cleanup_normaliz_files,
    get_all_unique_domains,
    load_covering_blocks,
    load_monomers,
    probe_k,
)
from paths import EXAMPLE_TBNS_DIR, RESULTS_DIR

PHASE2_JSON = os.path.join(RESULTS_DIR, "experiments", "phase2_covering.json")
K_FIXED     = 25


def _monomer_path(family, size):
    if family == "cascade":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{size}.txt")
    if family == "binary":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_binary_tree_d{size}.txt")
    if family == "dna":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_dna_tbn_depth{size}.txt")
    raise ValueError(family)


def main():
    with open(PHASE2_JSON) as f:
        data = json.load(f)

    targets = [r for r in data["results"] if r.get("best_k") is None]
    print(f"Refill targets: {len(targets)} cells")
    for r in targets:
        print(f"  {r['family']:8s} size={r['size']} t={r['t']}")

    if not targets:
        print("Nothing to refill.")
        return

    class _NullLog:
        def write(self, *a): pass
        def flush(self): pass
    log = _NullLog()

    for r in targets:
        family, size, t = r["family"], r["size"], r["t"]
        path = _monomer_path(family, size)
        if not os.path.exists(path):
            print(f"  [{family} size={size} t={t}] missing monomer file -- skip")
            continue
        all_monomers = load_monomers(path)
        all_domains  = get_all_unique_domains(all_monomers)
        n            = len(all_monomers)

        fallback_dp = (family == "binary" and size == 4)

        print(f"\n  [{family} size={size} t={t}]  |M|={n}  "
              f"fallback_dp={fallback_dp}  k={K_FIXED}")
        cleanup_normaliz_files()

        t0 = time.time()
        try:
            blocks = load_covering_blocks(n, K_FIXED, t, fallback_dp=fallback_dp)
        except Exception as e:
            print(f"    load_covering_blocks failed: {type(e).__name__}: {e}")
            r["best_k"] = None
            r["best_projected"] = None
            r["full"] = None
            r["full_run_skipped_reason"] = f"refill_failed: {e}"
            continue
        load_secs = time.time() - t0
        print(f"    covering loaded in {load_secs:.1f}s  ({len(blocks)} blocks)")

        t1 = time.time()
        projected, probe_times, num_blocks = probe_k(
            K_FIXED, t, blocks, n, all_monomers, "monomer", all_domains, n,
            best_projected=None, log=log,
        )
        probe_secs = time.time() - t1
        print(f"    probed {len(probe_times)}/{num_blocks} blocks in "
              f"{probe_secs:.1f}s -- projected={projected:.1f}s")

        r["best_k"] = K_FIXED
        r["best_projected"] = projected
        r["probe_summary"] = {str(K_FIXED): {
            "projected_total": projected,
            "num_blocks":      num_blocks,
            "probe_count":     len(probe_times),
        }}
        r["full"] = None
        r["full_run_skipped_reason"] = "refilled_probe_only"
        r["timestamp_refilled"] = datetime.now().isoformat()

        with open(PHASE2_JSON, "w") as f:
            json.dump(data, f, indent=2, default=str)

    print(f"\nDone -- phase2_covering.json updated.")


if __name__ == "__main__":
    main()
