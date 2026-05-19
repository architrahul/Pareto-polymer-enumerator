"""
leakage_experiment_inputs.py — EXPERIMENT 2.

Compare actual leakage across cascade-N variants where the first K inputs
are removed, for K = 1, ..., N. For each variant we compute the full
Hilbert basis P* (no covering), run COFFEE under the leakage scenario,
and compare per-polymer equilibrium concentrations against the
user-specified "expected" polymer set (see expected_polymers.py).

Outputs:
  results/benchmarks/04_leakage/removed_inputs/
  reusable shared cache: results/common/hilbert_basis/ and results/common/coffee/
    K{K}/
      coffee_output.txt
      input.ocx / input.con (in coffee_input/)
      polymers_sorted.csv          (top-N polymers + concentrations)
    summary.json                   (aggregate leakage per K)
    figure_leakage_vs_K.png        (bar chart)
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime


# Keep matplotlib cache inside results/ so scripts run cleanly on locked-down machines.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results", ".matplotlib"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coffee_parser import (
    assign_domain_energies,
    generate_coffee_inputs,
    modify_output,
    parse_monomers,
    read_polymers,
)
from hilbert_pipeline import (
    cleanup_normaliz_files,
    get_all_unique_domains,
    run_normaliz_on_subset,
    save_polymer_vectors,
)
from paths import (
    EXAMPLE_TBNS_DIR,
    LEAKAGE_HB_DIR,
    LEAKAGE_COFFEE_DIR,
    LEAKAGE_VARY_INPUTS_DIR,
    REPO_ROOT,
)

from expected_polymers import (
    DEFAULT_INITIAL_CONC,
    build_expected_polymers,
    generate_incomplete_cascade,
)

COFFEE_CLI = os.path.join(REPO_ROOT, "coffee", "crates", "coffee-cli",
                          "target", "release", "coffee-cli")
K_VALUES = None  # by default, use 1..n for the requested cascade depth

SIGNIFICANCE_CUTOFF = 1e-9   # 1 nM (1 µM initial monomer conc)


# ---------------------------------------------------------------------------
# Full-HB computation, with caching
# ---------------------------------------------------------------------------

def _compute_full_hb(monomers_path: str, system_label: str) -> str:
    """Run one Normaliz call on the full system; cache the projected basis at
    results/common/hilbert_basis/hilbert_full_p_star_{system_label}.txt.
    Returns the path (or raises if Normaliz returned no vectors)."""
    out = os.path.join(LEAKAGE_HB_DIR, f"hilbert_full_p_star_{system_label}.txt")
    if os.path.exists(out):
        print(f"  [{system_label}] Full P* cached -> {os.path.basename(out)}")
        return out

    with open(monomers_path) as f:
        all_monomers = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    n = len(all_monomers)
    print(f"  [{system_label}] running Full HB Normaliz on |M|={n} ...")
    cleanup_normaliz_files()
    t0 = time.time()
    elapsed, raw = run_normaliz_on_subset(all_monomers)
    print(f"    Normaliz finished in {elapsed:.1f}s, {len(raw)} raw vectors")
    if not raw:
        raise RuntimeError(f"Normaliz returned no vectors for {system_label}")

    # Project to monomer space (drop augmented unit-monomer columns), de-dup
    vecs = {tuple(v[:n]) for v in raw}
    vecs.discard(tuple([0] * n))
    os.makedirs(LEAKAGE_HB_DIR, exist_ok=True)
    save_polymer_vectors(
        vecs, out, n_monomers=n,
        comment=f"Full P* (Full HB), system={system_label}",
    )
    print(f"    saved {len(vecs)} vectors -> {os.path.basename(out)}")
    return out


# ---------------------------------------------------------------------------
# COFFEE + comparison helpers (shared with experiment 1)
# ---------------------------------------------------------------------------

def _run_coffee(ocx_path, con_path, out_path):
    assert os.path.isfile(COFFEE_CLI), f"coffee-cli not built at {COFFEE_CLI}"
    t0 = time.time()
    subprocess.run([COFFEE_CLI, ocx_path, con_path, "-o", out_path],
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False)
    return time.time() - t0


def _read_coffee_concs(out_path, n_polymers):
    with open(out_path) as f:
        for line in f:
            parts = line.strip().split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) == n_polymers:
                return np.array(vals)
    raise RuntimeError(f"could not parse {n_polymers} concentrations from {out_path}")


def _compare(polymers, concs, expected_polys):
    exp_by_vec = {tuple(ep.vector): ep for ep in expected_polys}
    rows = []
    for p, c in zip(polymers, concs):
        key = tuple(p)
        ep  = exp_by_vec.get(key)
        rows.append({
            "vector":         key,
            "conc":           float(c),
            "is_expected":    ep is not None,
            "expected_label": ep.label if ep else None,
            "expected_conc":  float(ep.concentration) if ep else 0.0,
            "deviation":      float(c) - (float(ep.concentration) if ep else 0.0),
        })
    return rows


def _aggregate(rows, expected_polys, cutoff=SIGNIFICANCE_CUTOFF):
    actual_for_exp = {ep.label: 0.0 for ep in expected_polys}
    unexpected = 0.0
    n_unexp_above = 0
    n_unexp_below = 0
    for r in rows:
        if r["is_expected"]:
            actual_for_exp[r["expected_label"]] = r["conc"]
        else:
            if r["conc"] >= cutoff:
                unexpected     += r["conc"]
                n_unexp_above  += 1
            elif r["conc"] > 0:
                n_unexp_below  += 1
    abs_dev = unexpected
    deficit = 0.0
    rec = miss = 0
    for ep in expected_polys:
        a_raw = actual_for_exp[ep.label]
        a = a_raw if a_raw >= cutoff else 0.0
        abs_dev += abs(a - ep.concentration)
        if a < ep.concentration:
            deficit += ep.concentration - a
        if a > 0:
            rec += 1
        else:
            miss += 1
    return {
        "cutoff":                 cutoff,
        "total_unexpected_conc":  unexpected,
        "total_expected_deficit": deficit,
        "total_abs_deviation":    abs_dev,
        "expected_recovered":     rec,
        "expected_missing":       miss,
        "n_unexpected_above":     n_unexp_above,
        "n_unexpected_below":     n_unexp_below,
        "expected_total":         len(expected_polys),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7.")
    parser.add_argument("--K-values", nargs="+", type=int, default=K_VALUES,
                        dest="k_values",
                        help="Inputs-removed counts. Default: 1..n.")
    parser.add_argument("--top", type=int, default=50,
                        help="Top-N polymers in the sorted dump. Default 50.")
    args = parser.parse_args()

    n        = args.cascade_n
    k_values = args.k_values or list(range(1, n + 1))
    out_root = os.path.join(LEAKAGE_VARY_INPUTS_DIR,
                            f"n{n}_systems_compare")
    os.makedirs(out_root, exist_ok=True)

    summary = {"generated": datetime.now().isoformat(), "by_K": {}}

    for K in k_values:
        sys_label = f"n{n}_incomplete{K}" if K > 1 else f"n{n}_incomplete"
        print(f"\n{'='*70}\n  System: cascade-{n}, first {K} input(s) removed "
              f"-> {sys_label}\n{'='*70}")

        # 1) Monomer file
        monomer_path = generate_incomplete_cascade(n, K)
        monomers      = parse_monomers(monomer_path)
        n_monomers    = len(monomers)
        domain_energy = assign_domain_energies(monomers, seed=42)
        # The first K inputs are removed from the monomer file entirely; the
        # surviving inputs x_{K+1}..x_{n+1} are kept at the initial 1 µM
        # concentration. We do NOT zero them in COFFEE.
        print(f"  |M|={n_monomers}  (all monomers at 1 µM)")

        # 2) Full HB
        try:
            hb_path = _compute_full_hb(monomer_path, sys_label)
        except Exception as e:
            print(f"  [{sys_label}] Full HB FAILED: {e}")
            continue

        # 3) Expected polymer set
        expected = build_expected_polymers(n, removed_inputs=K,
                                           initial_conc=DEFAULT_INITIAL_CONC)
        print(f"  expected polymer set: {len(expected)} polymers")

        # 4) COFFEE on full P*. Keep this experiment's COFFEE outputs separate
        # from the paper Figure 5/6 COFFEE outputs: those zero all input-signal
        # monomers, whereas here the removed inputs are absent from the monomer
        # file and every surviving monomer stays at 1 µM.
        sys_out_dir = os.path.join(out_root, f"K{K}")
        os.makedirs(sys_out_dir, exist_ok=True)
        shared_cof = os.path.join(LEAKAGE_COFFEE_DIR, "vary_removed_inputs",
                                  sys_label, "full_pstar",
                                  "coffee_output.txt")
        if os.path.exists(shared_cof):
            cof_dir = os.path.dirname(shared_cof)
            ocx = os.path.join(cof_dir, "input.ocx")
            con = os.path.join(cof_dir, "input.con")
            cof = shared_cof
            print(f"    COFFEE cached at {cof}")
        else:
            cof_dir = os.path.join(LEAKAGE_COFFEE_DIR, "vary_removed_inputs",
                                   sys_label, "full_pstar")
            os.makedirs(cof_dir, exist_ok=True)
            ocx = os.path.join(cof_dir, "input.ocx")
            con = os.path.join(cof_dir, "input.con")
            cof = os.path.join(cof_dir, "coffee_output.txt")
            generate_coffee_inputs(monomers=monomers, polymer_file=hb_path,
                                   out_dir=cof_dir,
                                   domain_energy=domain_energy,
                                   zero_indices=None,  # surviving inputs at 1 µM
                                   label=f"full_pstar_K{K}")
            elapsed = _run_coffee(ocx, con, cof)
            print(f"    coffee-cli on full P* in {elapsed:.1f}s")

        csv_out = modify_output(cof, hb_path)
        print(f"    wrote {os.path.basename(csv_out)}")

        # 5) Per-polymer comparison
        polymers = read_polymers(hb_path)
        concs    = _read_coffee_concs(cof, len(polymers))
        rows     = _compare(polymers, concs, expected)
        rows.sort(key=lambda r: -r["conc"])

        # 6) Sorted polymer dump — cutoff-filtered (instead of fixed top-N).
        # CSV: only above-cutoff polymers + all expected (so suppressed
        # expected ones are still visible to the user).
        csv_path = os.path.join(sys_out_dir, "polymers_sorted.csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "concentration_M", "above_cutoff",
                        "is_expected", "expected_label",
                        "expected_conc_M", "deviation_M",
                        "vector_nonzero_indices", "full_polymer_vector"])
            rank = 0
            for r in rows:
                above = r["conc"] >= SIGNIFICANCE_CUTOFF
                if not above and not r["is_expected"]:
                    continue
                nz = " ".join(str(j) for j, v in enumerate(r["vector"]) if v)
                full_vec = " ".join(str(v) for v in r["vector"])
                w.writerow([rank, f"{r['conc']:.6e}", above,
                            r["is_expected"],
                            r["expected_label"] or "",
                            f"{r['expected_conc']:.6e}",
                            f"{r['deviation']:.6e}", nz, full_vec])
                rank += 1

        # Human-readable: every polymer above the cutoff (any count)
        top_path = os.path.join(
            sys_out_dir,
            f"polymers_above_{SIGNIFICANCE_CUTOFF:.0e}.txt",
        )
        rows_above = [r for r in rows if r["conc"] >= SIGNIFICANCE_CUTOFF]
        with open(top_path, "w") as f:
            f.write(f"# Polymers with equilibrium conc >= "
                    f"{SIGNIFICANCE_CUTOFF:g} M  "
                    f"(cascade-{n}, {K} input(s) removed)\n")
            f.write(f"# total polymers in P*: {len(rows)}\n")
            f.write(f"# above cutoff: {len(rows_above)}\n")
            f.write(f"# expected polymer set size: {len(expected)}\n\n")
            for i, r in enumerate(rows_above):
                tag = f"EXPECTED:{r['expected_label']}" if r["is_expected"] else "unexpected"
                monomer_strs = [
                    " ".join(monomers[j]) for j, v in enumerate(r["vector"]) if v
                ]
                f.write(f"{i:>4}  c={r['conc']:.3e} M  [{tag}]  :: "
                        f"{' + '.join(monomer_strs)}\n")
                f.write(f"      vector: {' '.join(str(v) for v in r['vector'])}\n")

        agg = _aggregate(rows, expected)
        summary["by_K"][K] = {
            "system_label":      sys_label,
            "n_monomers":        n_monomers,
            "n_expected":        len(expected),
            "n_full_pstar":      len(polymers),
            "aggregate":         agg,
        }
        print(f"  recovered={agg['expected_recovered']}/{len(expected)}  "
              f"unexpected_conc={agg['total_unexpected_conc']:.3e}  "
              f"abs_dev={agg['total_abs_deviation']:.3e}")

    with open(os.path.join(out_root, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    _plot_summary(summary, out_root)


def _plot_summary(summary, out_dir):
    Ks  = sorted(int(k) for k in summary["by_K"])
    if not Ks:
        return
    metrics = ("total_unexpected_conc",
               "total_expected_deficit",
               "total_abs_deviation")

    fig, ax = plt.subplots(figsize=(8, 5))
    x     = np.arange(len(Ks))
    bar_w = 0.27
    for i, m in enumerate(metrics):
        heights = [summary["by_K"][k]["aggregate"][m] for k in Ks]
        off = (i - len(metrics) / 2 + 0.5) * bar_w
        ax.bar(x + off, heights, bar_w, label=m,
               edgecolor="black", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels([f"K={k}" for k in Ks])
    ax.set_yscale("log")
    ax.set_ylabel("concentration (M) — log")
    ax.set_xlabel("number of first inputs removed")
    ax.set_title("Actual leakage (full P*) across cascade systems")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, which="both", axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    out = os.path.join(out_dir, "figure_leakage_vs_K.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
