"""
leakage_experiment_t.py — EXPERIMENT 1.

Worst-case leakage analysis: cascade-N with the first input (x1) removed.
For each t in {3..8} we compute the covering-design enumeration P̂_{k=25,t},
plus the full Pareto-optimal P* (already computed) and the user-defined
"expected" polymer set (see expected_polymers.py). COFFEE is run on each
under the leakage scenario (x1 = 0 µM, others = 1 µM). We then compare:

    actual_leakage  = full P*       equilibrium  vs  expected
    analysed_leakage(t) = P̂_{k=25,t} equilibrium  vs  expected

across t, and write:
    results/benchmarks/04_leakage/vary_t/
  reusable shared cache: results/common/hilbert_basis/ and results/common/coffee/
        figure_leakage_vs_t.png           summary plot
        summary.json                      aggregate metrics
        polymer_compare_full_pstar.csv    per-polymer table, full P*
        polymer_compare_k25_t{T}.csv      per-polymer table, P̂_{k=25,t}
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
from paths import (
    EXAMPLE_TBNS_DIR,
    LEAKAGE_HB_DIR,
    LEAKAGE_COFFEE_DIR,
    LEAKAGE_VARY_T_DIR,
    REPO_ROOT,
)

from expected_polymers import (
    DEFAULT_INITIAL_CONC,
    build_expected_polymers,
    generate_incomplete_cascade,
)
from leakage_analysis import leakage_analysis

COFFEE_CLI = os.path.join(REPO_ROOT, "coffee", "crates", "coffee-cli",
                          "target", "release", "coffee-cli")
T_VALUES   = [3, 4, 5, 6, 7]   # t=8 skipped — covering compute takes ~25 min
K_FIXED    = 25

# No concentration cutoff for leakage totals: every expected and unexpected
# candidate polymer contributes to the concentration-error metrics.
SIGNIFICANCE_CUTOFF = 0.0


# ---------------------------------------------------------------------------
# COFFEE helpers
# ---------------------------------------------------------------------------

def _run_coffee(ocx_path, con_path, out_path):
    assert os.path.isfile(COFFEE_CLI), f"coffee-cli not built at {COFFEE_CLI}"
    t0 = time.time()
    subprocess.run([COFFEE_CLI, ocx_path, con_path, "-o", out_path],
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False)
    return time.time() - t0


def _read_coffee_concs(out_path, n_polymers):
    """First line of `out_path` with exactly n_polymers floats."""
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


def _save_expected_polymer_file(expected_polys, out_path):
    """Write expected polymer vectors as a polymer-file (one vector per line)."""
    with open(out_path, "w") as f:
        f.write(f"# Expected polymer set under leakage  "
                f"({len(expected_polys)} polymers, init conc 1 µM)\n")
        for ep in expected_polys:
            f.write(" ".join(str(x) for x in ep.vector) + "\n")


# ---------------------------------------------------------------------------
# Per-polymer leakage comparison
# ---------------------------------------------------------------------------

def _compare(polymers, concs, expected_polys, initial_conc):
    """For each polymer in `polymers`, decide if it matches an expected one
    (by vector equality) and compute its concentration relative to expected.

    Returns a list of dicts with keys: rank, vector, conc, is_expected,
    expected_label (or None), expected_conc (0 if not expected),
    deviation = conc - expected_conc.
    """
    exp_by_vec = {tuple(ep.vector): ep for ep in expected_polys}
    rows = []
    for p, c in zip(polymers, concs):
        key = tuple(p)
        ep  = exp_by_vec.get(key)
        rows.append({
            "vector":        key,
            "conc":          float(c),
            "is_expected":   ep is not None,
            "expected_label": ep.label if ep else None,
            "expected_conc": float(ep.concentration) if ep else 0.0,
            "deviation":     float(c) - (float(ep.concentration) if ep else 0.0),
        })
    return rows


def _aggregate_metrics(rows, expected_polys, cutoff=SIGNIFICANCE_CUTOFF):
    """Aggregate leakage metrics across a single polymer set's equilibrium.

    Polymers with actual concentration >= `cutoff` contribute. In the current
    leakage experiment cutoff=0, so all candidate polymers contribute.

    Returns:
      cutoff:                 the threshold used
      total_unexpected_conc:  sum of conc over UNEXPECTED polymers with conc >= cutoff
      total_expected_deficit: sum over expected polymers of max(0, expected_conc - max(actual, 0))
                              — i.e. for expected polymers below cutoff we use actual=0
      total_abs_deviation:    sum of |actual - expected| across both groups
                              (unexpected below cutoff contributes 0; expected below
                               cutoff contributes full expected_conc)
      expected_recovered:     count of expected polymers with actual >= cutoff
      expected_missing:       count of expected polymers with actual < cutoff
      n_unexpected_above:     count of unexpected polymers with conc >= cutoff
      n_unexpected_below:     count of unexpected polymers with 0 < conc < cutoff (informational)
    """
    actual_for_expected = {ep.label: 0.0 for ep in expected_polys}
    unexpected_conc = 0.0
    n_unexp_above = 0
    n_unexp_below = 0
    for r in rows:
        if r["is_expected"]:
            actual_for_expected[r["expected_label"]] = r["conc"]
        else:
            if r["conc"] >= cutoff:
                unexpected_conc += r["conc"]
                n_unexp_above   += 1
            elif r["conc"] > 0:
                n_unexp_below   += 1

    abs_dev = unexpected_conc
    deficit = 0.0
    recovered = missing = 0
    for ep in expected_polys:
        a_raw = actual_for_expected[ep.label]
        # Treat sub-cutoff actuals as 0 for deficit/recovery purposes
        a = a_raw if a_raw >= cutoff else 0.0
        abs_dev += abs(a - ep.concentration)
        if a < ep.concentration:
            deficit += ep.concentration - a
        if a > 0:
            recovered += 1
        else:
            missing += 1

    return {
        "cutoff":                 cutoff,
        "total_unexpected_conc":  unexpected_conc,
        "total_expected_deficit": deficit,
        "total_abs_deviation":    abs_dev,
        "expected_recovered":     recovered,
        "expected_missing":       missing,
        "n_unexpected_above":     n_unexp_above,
        "n_unexpected_below":     n_unexp_below,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7.")
    parser.add_argument("--logs", action="store_true",
                        help="Write verbose covering-enumeration logs. Default: no verbose log files.")
    args = parser.parse_args()

    n = args.cascade_n
    label = f"n{n}_incomplete"
    out_dir = os.path.join(LEAKAGE_VARY_T_DIR, label)
    os.makedirs(out_dir, exist_ok=True)

    monomer_file = generate_incomplete_cascade(n, removed_inputs=1)
    monomers      = parse_monomers(monomer_file)
    n_monomers    = len(monomers)
    domain_energy = assign_domain_energies(monomers, seed=42)
    # NOTE: the "removed" input (x1) is gone from the monomer file entirely.
    # All remaining monomers — including the surviving free inputs x2..x8 —
    # are at the initial 1 µM concentration. No COFFEE-side zeroing here.
    print(f"System: cascade {label}, |M|={n_monomers}  (all monomers at 1 µM)")

    # 1) Build & persist expected polymer set
    expected = build_expected_polymers(n, removed_inputs=1,
                                       initial_conc=DEFAULT_INITIAL_CONC)
    print(f"Expected polymer set: {len(expected)} polymers")
    expected_path = os.path.join(out_dir, "expected_polymers.txt")
    _save_expected_polymer_file(expected, expected_path)

    # 2) Ensure P̂_{k=25, t} files exist for all t (run covering for missing t)
    print()
    for t in T_VALUES:
        out_txt = os.path.join(
            LEAKAGE_HB_DIR,
            f"hilbert_k{K_FIXED}_t{t}_monomer_{label}.txt",
        )
        if os.path.exists(out_txt):
            print(f"  [covering t={t}] exists  -> {os.path.basename(out_txt)}")
            continue
        print(f"  [covering t={t}] missing -> running leakage_analysis ...")
        leakage_analysis(cascade_n=n, t=t, k=K_FIXED, only="incomplete", logs=args.logs)

    # 3) Discover all polymer files we want to compare.
    # NOTE: the expected polymer set is NOT run through COFFEE — it's the
    # reference by construction (each polymer at the initial monomer conc).
    # COFFEE requires #polymers >= #monomers anyway; the 28 expected
    # polymers vs 56 monomers would fail the optimizer.
    polymer_sources = {}
    polymer_sources["full_pstar"] = os.path.join(
        LEAKAGE_HB_DIR, f"hilbert_full_p_star_{label}.txt"
    )
    for t in T_VALUES:
        polymer_sources[f"k{K_FIXED}_t{t}"] = os.path.join(
            LEAKAGE_HB_DIR, f"hilbert_k{K_FIXED}_t{t}_monomer_{label}.txt"
        )

    # 4) Build COFFEE inputs + run COFFEE for each set. Reuse the canonical
    # shared cache results/common/coffee/{system}/{polymer_set}/ so experiment
    # 3 and experiment 4 do not recompute identical COFFEE analyses.
    shared_coffee_base = os.path.join(LEAKAGE_COFFEE_DIR, label)
    os.makedirs(shared_coffee_base, exist_ok=True)
    coffee_results = {}

    for tag, p_path in polymer_sources.items():
        if not os.path.exists(p_path):
            print(f"  [{tag}] polymer file missing: {p_path} - SKIP")
            continue
        d = os.path.join(shared_coffee_base, tag)
        ocx = os.path.join(d, "input.ocx")
        con = os.path.join(d, "input.con")
        cof = os.path.join(d, "coffee_output.txt")

        if os.path.exists(cof):
            print(f"  [{tag}] COFFEE output cached")
        else:
            print(f"  [{tag}] building COFFEE inputs + running coffee-cli ...")
            generate_coffee_inputs(
                monomers=monomers,
                polymer_file=p_path,
                out_dir=d,
                domain_energy=domain_energy,
                zero_indices=None,    # nothing zeroed — all surviving monomers at 1 µM
                label=tag,
            )
            elapsed = _run_coffee(ocx, con, cof)
            print(f"    coffee-cli in {elapsed:.1f}s")

        csv_out = modify_output(cof, p_path)
        print(f"    wrote {os.path.basename(csv_out)}")

        polymers = read_polymers(p_path)
        concs    = _read_coffee_concs(cof, len(polymers))
        coffee_results[tag] = (polymers, concs)

    # 5) Per-polymer leakage tables. Compared against the constructed expected
    # set (no COFFEE needed for the expected reference). Leakage totals use no
    # concentration cutoff: every expected and unexpected candidate contributes.
    summary = {"generated": datetime.now().isoformat(),
               "n_expected_polymers": len(expected),
               "significance_cutoff_M": SIGNIFICANCE_CUTOFF,
               "by_set": {}}
    for tag, (polymers, concs) in coffee_results.items():
        rows = _compare(polymers, concs, expected, DEFAULT_INITIAL_CONC)
        rows.sort(key=lambda r: -r["conc"])

        # CSV: include all candidate polymers (cutoff=0), plus expected labels.
        csv_path = os.path.join(out_dir, f"polymer_compare_{tag}.csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "concentration_M", "above_cutoff",
                        "is_expected", "expected_label",
                        "expected_conc_M", "deviation_M",
                        "vector_nonzero_indices", "full_polymer_vector"])
            rank = 0
            for r in rows:
                above = r["conc"] >= SIGNIFICANCE_CUTOFF
                # With cutoff=0, all candidate polymers are included.
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

        agg = _aggregate_metrics(rows, expected)
        summary["by_set"][tag] = agg
        print(f"  [{tag}] |P|={len(polymers)}  "
              f"recovered={agg['expected_recovered']}/{len(expected)} (>=1nM)  "
              f"unexpected_conc(>=1nM)={agg['total_unexpected_conc']:.3e}  "
              f"abs_dev={agg['total_abs_deviation']:.3e}  "
              f"n_unexp_above={agg['n_unexpected_above']}")

    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # 6) Summary plot
    _plot_summary(summary, out_dir, label)


def _plot_summary(summary, out_dir, label):
    sets = list(summary["by_set"].keys())
    metrics = ("total_unexpected_conc",
               "total_expected_deficit",
               "total_abs_deviation")

    # Reorder so full_pstar is leftmost, then k=25,t in ascending t
    def _key(tag):
        if tag == "full_pstar":
            return (0, 0)
        if tag.startswith("k25_t"):
            return (1, int(tag.split("t")[-1]))
        return (2, tag)
    sets.sort(key=_key)

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(sets))
    bar_w = 0.27

    for i, m in enumerate(metrics):
        heights = [summary["by_set"][s][m] for s in sets]
        off     = (i - len(metrics) / 2 + 0.5) * bar_w
        ax.bar(x + off, heights, bar_w, label=m, edgecolor="black", linewidth=0.3)

    full_abs = summary["by_set"].get("full_pstar", {}).get("total_abs_deviation")
    if full_abs is not None:
        ax.axhline(full_abs, color="black", linestyle="--", linewidth=0.9,
                   alpha=0.7, label="Full P* abs. deviation")

    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", " ") for s in sets], rotation=0)
    ax.set_ylim(0, max(summary["by_set"][s][m] for s in sets for m in metrics) * 1.12)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_ylabel("concentration (M)")
    ax.set_title(f"Leakage analysis — cascade {label}, leakage scenario")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    out = os.path.join(out_dir, "figure_leakage_vs_t.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
