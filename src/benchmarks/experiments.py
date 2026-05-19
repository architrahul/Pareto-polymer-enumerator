"""
experiments.py — Orchestrator for the paper's Figure 2/3/4 benchmarks.

Phase 1 (Figure 2): Probe-estimated total runtime as a function of block size k,
  t=5, for damien (|M|=21), cascade n=7 (|M|=57), cascade n=8 (|M|=65).
  Saves per-k projected-total to JSON.

Phase 2 (Figure 3 + Figure 4 covering bars): Probe + full enumeration on the
  best k, for cascade m=5..9, binary tree d=3..4, dna cascade m=4..7,
  t in {3..8}. Saves per-cell timing/vector counts as JSON.

Phase 3 (Figure 3 + Figure 4 "Full HB" bars): Full Hilbert basis on each
  (system, size) once (single Normaliz call on full system). The 3-hour
  per-Normaliz timeout enforced inside run_normaliz_on_subset will produce
  "hrs+" truncation bars matching the paper.

Phase 5 (Figure 5/6): leakage analysis on cascade n=7 and n=8 — already run
  by run_all_benchmarks.sh; this script does not repeat it.
"""

import argparse
import json
import os
import sys
import time
import traceback
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # parent src/

from paths import EXAMPLE_TBNS_DIR, RESULTS_DIR, LOGS_DIR
from hilbert_pipeline import (
    K_MAX,
    NORMALIZ_TIMEOUT_SECONDS,
    cleanup_normaliz_files,
    full_run_k,
    get_all_unique_domains,
    load_covering_blocks,
    load_monomers,
    probe_k,
    run_normaliz_on_subset,
    save_polymer_vectors,
    start_input_listener,
)

EXP_DIR    = os.path.join(RESULTS_DIR, "experiments")
HB_OUT_DIR = os.path.join(EXP_DIR, "hilbert_basis")
os.makedirs(EXP_DIR, exist_ok=True)
os.makedirs(HB_OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# System registry — paper-aligned
# ---------------------------------------------------------------------------

CASCADE_SIZES = [5, 6, 7, 8, 9]
BINARY_SIZES  = [3, 4]
DNA_SIZES     = [4, 5, 6, 7]

T_VALUES = [3, 4, 5, 6, 7, 8]

# Figure 2 systems: paper Section 6.2 uses damien (scaffolded DNA, |M|=21) and
# the 8-module linear cascade. Because the paper's "n=8 with |M|=57" is the
# typo, we run both n=7 and n=8 here.
FIGURE2_SYSTEMS = [
    ("damien_n10", "monomers_damien_n10.txt"),
    ("cascade_n7", "monomers_cascade_n7.txt"),
    ("cascade_n8", "monomers_cascade_n8.txt"),
]

FIGURE2_T = 5

# Phase 2 only does a full enumeration on the best k if the *projected* runtime
# for that k is below this cap. Larger cells just record the probe-derived
# projection (which is informative enough for Figure 3/4 bars).
PHASE2_FULL_RUN_CAP_S = 30 * 60


def _path(family: str, size: int) -> str:
    if family == "cascade":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{size}.txt")
    if family == "binary":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_binary_tree_d{size}.txt")
    if family == "dna":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_dna_tbn_depth{size}.txt")
    raise ValueError(family)


def _save_json(path, payload):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# Phase 1 — Figure 2 probe curves
# ---------------------------------------------------------------------------

def phase1_probe_curves(log):
    """Probe every k in [t+1, K_MAX] (plus k=n) and record projected totals.

    No pruning — we want the full curve, not the argmin.
    """
    print("\n" + "=" * 70)
    print("PHASE 1 — Figure 2 probe curves (t=5)")
    print("=" * 70)

    for label, fname in FIGURE2_SYSTEMS:
        monomer_file = os.path.join(EXAMPLE_TBNS_DIR, fname)
        all_monomers = load_monomers(monomer_file)
        n            = len(all_monomers)
        all_domains  = get_all_unique_domains(all_monomers)
        t            = FIGURE2_T

        out_path = os.path.join(EXP_DIR, f"figure2_{label}.json")
        print(f"\n  [{label}] |M|={n}, t={t}  →  {out_path}")
        log.write(f"\n[phase1] {label}: |M|={n}, t={t}\n"); log.flush()

        k_values = list(range(t + 1, min(K_MAX, n) + 1))
        if n not in k_values:
            k_values.append(n)

        payload = {
            "label":         label,
            "monomer_file":  monomer_file,
            "n_monomers":    n,
            "n_domains":     len(all_domains),
            "t":             t,
            "k_values":      k_values,
            "results":       {},
            "started":       datetime.now().isoformat(),
        }
        _save_json(out_path, payload)

        for k in k_values:
            cleanup_normaliz_files()
            t0 = time.time()
            try:
                blocks = load_covering_blocks(n, k, t, fallback_dp=False)
            except Exception as e:
                payload["results"][str(k)] = {"error": f"load_covering_blocks failed: {e}"}
                _save_json(out_path, payload); continue

            projected, probe_times, num_blocks = probe_k(
                k, t, blocks, n, all_monomers, "monomer", all_domains, n,
                best_projected=None, log=log,
            )

            payload["results"][str(k)] = {
                "k":               k,
                "num_blocks":      num_blocks,
                "projected_total": projected,
                "probe_count":     len(probe_times),
                "probe_total":     sum(probe_times),
                "wall":            time.time() - t0,
            }
            _save_json(out_path, payload)
            print(f"    k={k:2d}  blocks={num_blocks:5d}  projected={projected}")

        payload["finished"] = datetime.now().isoformat()
        _save_json(out_path, payload)


# ---------------------------------------------------------------------------
# Phase 2 — covering bars (Figures 3 and 4)
# ---------------------------------------------------------------------------

def phase2_covering(log, families):
    """For each (family, size, t), probe+full enumerate on the best k."""
    print("\n" + "=" * 70)
    print("PHASE 2 — covering sweep (Figures 3 & 4 covering bars)")
    print("=" * 70)

    cells = []
    if "cascade" in families:
        cells += [("cascade", s, t) for s in CASCADE_SIZES for t in T_VALUES]
    if "binary" in families:
        cells += [("binary",  s, t) for s in BINARY_SIZES  for t in T_VALUES]
    if "dna" in families:
        cells += [("dna",     s, t) for s in DNA_SIZES     for t in T_VALUES]

    results_path = os.path.join(EXP_DIR, "phase2_covering.json")
    results = []
    if os.path.exists(results_path):
        try:
            results = json.load(open(results_path))["results"]
        except Exception:
            results = []
    done_keys = {(r["family"], r["size"], r["t"]) for r in results}

    for idx, (family, size, t) in enumerate(cells, 1):
        if (family, size, t) in done_keys:
            print(f"  [{idx}/{len(cells)}] {family} size={size} t={t}  SKIP (cached)")
            continue

        path = _path(family, size)
        if not os.path.exists(path):
            log.write(f"[phase2] missing monomer file: {path}\n"); log.flush()
            continue

        all_monomers = load_monomers(path)
        n            = len(all_monomers)
        all_domains  = get_all_unique_domains(all_monomers)

        if t + 1 > min(K_MAX, n):
            log.write(f"[phase2] {family} size={size} t={t}: SKIP (no valid k)\n"); log.flush()
            continue

        print(f"\n  [{idx}/{len(cells)}] {family} size={size} t={t}  (|M|={n})")
        log.write(f"\n[phase2] {family} size={size} t={t} |M|={n}\n"); log.flush()
        cell_start = time.time()

        # ---- probe sweep ----
        best_k, best_projected, best_blocks = None, None, None
        probe_summary = {}
        for k in range(t + 1, min(K_MAX, n) + 1):
            cleanup_normaliz_files()
            try:
                blocks = load_covering_blocks(n, k, t, fallback_dp=False)
            except Exception as e:
                probe_summary[k] = {"error": str(e)}
                continue
            projected, probe_times, num_blocks = probe_k(
                k, t, blocks, n, all_monomers, "monomer", all_domains, n,
                best_projected=best_projected, log=log,
            )
            probe_summary[k] = {
                "projected_total": projected,
                "num_blocks":      num_blocks,
                "probe_count":     len(probe_times),
            }
            if projected is not None and (best_projected is None or projected < best_projected):
                best_projected, best_k, best_blocks = projected, k, blocks

        # ---- full run on best k (only if projected runtime <= cap) ----
        full = None
        full_run_skipped_reason = None
        if best_k is None:
            full_run_skipped_reason = "no_valid_k"
        elif best_projected > PHASE2_FULL_RUN_CAP_S:
            full_run_skipped_reason = f"projected_over_cap_{PHASE2_FULL_RUN_CAP_S}s"
            msg = (f"    best k={best_k}  projected={best_projected:.1f}s "
                   f"> cap={PHASE2_FULL_RUN_CAP_S}s — skipping full enumeration")
            print(msg)
            log.write(f"[phase2] {msg.strip()}\n"); log.flush()
        else:
            print(f"    best k={best_k}  projected={best_projected:.3f}s — running full enumeration")
            log.write(f"[phase2] best_k={best_k} projected={best_projected:.3f}s\n"); log.flush()
            cleanup_normaliz_files()
            full_result, _ = full_run_k(
                best_k, best_blocks, n, all_monomers, "monomer", all_domains, n, log
            )
            if full_result is not None:
                hb_path = os.path.join(
                    HB_OUT_DIR,
                    f"phase2_{family}_size{size}_t{t}_k{best_k}.txt"
                )
                save_polymer_vectors(
                    full_result["vectors"], hb_path,
                    n_monomers=n,
                    comment=f"phase2 covering, family={family} size={size} t={t} k={best_k}",
                )
                full = {
                    "k":                 best_k,
                    "wall":              full_result["total_wall_time"],
                    "normaliz_time":     full_result["total_normaliz_time"],
                    "overhead":          full_result["overhead_time"],
                    "unique_vectors":    full_result["unique_vectors"],
                    "hb_path":           hb_path,
                }

        results.append({
            "family":          family,
            "size":            size,
            "t":               t,
            "n_monomers":      n,
            "best_k":          best_k,
            "best_projected":  best_projected,
            "probe_summary":   probe_summary,
            "full":            full,
            "full_run_skipped_reason": full_run_skipped_reason,
            "cell_wall":       time.time() - cell_start,
            "timestamp":       datetime.now().isoformat(),
        })
        _save_json(results_path, {"generated": datetime.now().isoformat(), "results": results})


# ---------------------------------------------------------------------------
# Phase 3 — Full HB bars (Figure 3 + Figure 4)
# ---------------------------------------------------------------------------

def phase3_full_hb(log, families):
    """For each (family, size), run Normaliz on the full system once.

    With NORMALIZ_TIMEOUT_SECONDS = 3h, this produces the paper's "hrs+"
    truncated bars when Full HB doesn't complete.
    """
    print("\n" + "=" * 70)
    print("PHASE 3 — Full HB bars  (timeout per run: "
          f"{NORMALIZ_TIMEOUT_SECONDS//3600}h)")
    print("=" * 70)

    targets = []
    if "cascade" in families: targets += [("cascade", s) for s in CASCADE_SIZES]
    if "binary"  in families: targets += [("binary",  s) for s in BINARY_SIZES]
    if "dna"     in families: targets += [("dna",     s) for s in DNA_SIZES]

    results_path = os.path.join(EXP_DIR, "phase3_full_hb.json")
    results = []
    if os.path.exists(results_path):
        try:
            results = json.load(open(results_path))["results"]
        except Exception:
            results = []
    done_keys = {(r["family"], r["size"]) for r in results}

    for idx, (family, size) in enumerate(targets, 1):
        if (family, size) in done_keys:
            print(f"  [{idx}/{len(targets)}] {family} size={size}  SKIP (cached)")
            continue

        path = _path(family, size)
        if not os.path.exists(path):
            log.write(f"[phase3] missing monomer file: {path}\n"); log.flush()
            continue

        all_monomers = load_monomers(path)
        n            = len(all_monomers)

        print(f"\n  [{idx}/{len(targets)}] {family} size={size}  |M|={n}  →  Full HB")
        log.write(f"\n[phase3] {family} size={size} |M|={n}\n"); log.flush()

        cleanup_normaliz_files()
        wall_t0 = time.time()
        elapsed, raw_vectors = run_normaliz_on_subset(all_monomers)
        wall = time.time() - wall_t0
        truncated = (elapsed >= NORMALIZ_TIMEOUT_SECONDS - 1)

        if not truncated and raw_vectors:
            hb_path = os.path.join(
                HB_OUT_DIR, f"phase3_full_hb_{family}_size{size}.txt"
            )
            # Project augmented-system vectors (|M| + 2|Σ| coords) down to
            # just the M coordinates (π map of Theorem 9). Drops the zero
            # vector and any unit-monomer duplicates.
            vector_tuples = {tuple(v[:n]) for v in raw_vectors}
            vector_tuples.discard(tuple([0] * n))
            save_polymer_vectors(
                vector_tuples, hb_path, n_monomers=n,
                comment=f"phase3 Full HB, family={family} size={size}",
            )
        else:
            hb_path = None

        results.append({
            "family":          family,
            "size":            size,
            "n_monomers":      n,
            "wall":            wall,
            "normaliz_time":   elapsed,
            "vectors_found":   len(raw_vectors),
            "truncated":       truncated,
            "hb_path":         hb_path,
            "timestamp":       datetime.now().isoformat(),
        })
        _save_json(results_path, {"generated": datetime.now().isoformat(), "results": results})

        status = "TRUNCATED" if truncated else f"{elapsed:.1f}s ({len(raw_vectors)} vectors)"
        msg = (f"    Full HB on {family} size={size} |M|={n}: "
               f"normaliz={elapsed:.1f}s wall={wall:.1f}s status={status}")
        print(msg)
        log.write(msg + "\n"); log.flush()


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--phases", nargs="+",
                        default=["1", "2", "3"],
                        help="Which phases to run. Default: all.")
    parser.add_argument("--families", nargs="+",
                        default=["cascade", "binary", "dna"],
                        help="Which families for Phases 2/3.")
    args = parser.parse_args()

    os.makedirs(LOGS_DIR, exist_ok=True)
    log_path = os.path.join(LOGS_DIR, f"experiments_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    print(f"Log: {log_path}")

    start_input_listener()
    with open(log_path, "w") as log:
        log.write(f"experiments.py started {datetime.now().isoformat()}\n")
        log.write(f"phases={args.phases} families={args.families}\n")
        log.flush()

        for phase in args.phases:
            try:
                if phase == "1":
                    phase1_probe_curves(log)
                elif phase == "2":
                    phase2_covering(log, args.families)
                elif phase == "3":
                    phase3_full_hb(log, args.families)
            except KeyboardInterrupt:
                print(f"\n[phase {phase}] interrupted; moving on")
                log.write(f"\n[phase {phase}] interrupted at {datetime.now().isoformat()}\n")
                log.flush()
            except Exception:
                tb = traceback.format_exc()
                print(f"\n[phase {phase}] EXCEPTION:\n{tb}")
                log.write(f"\n[phase {phase}] EXCEPTION:\n{tb}\n"); log.flush()

        log.write(f"\nexperiments.py finished {datetime.now().isoformat()}\n")

    cleanup_normaliz_files()


if __name__ == "__main__":
    main()
