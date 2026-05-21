"""
export_csv.py — Export all numeric experiment results as tidy CSV files
into results/csv/, so external plotting tools can re-make the figures.

Reads from results/experiments/ and results/leakage/analysis/, writes:

  figure2_probe_curves.csv         (system, k) probe data — Figure 2
  figure3_4_covering_bars.csv      (family, size, t) — Figures 3/4 covering
  figure3_4_full_hb.csv            (family, size) — Figures 3/4 Full HB
  leakage_vs_t.csv                 Exp 1 aggregates per polymer set
  leakage_vs_t_per_polymer.csv     Exp 1 long-form per-polymer table
  leakage_vs_K.csv                 Exp 2 aggregates per K
  leakage_vs_K_per_polymer.csv     Exp 2 long-form per-polymer table

Every file has a CSV header on the first line.
"""

import csv
import json
import os
import sys
from glob import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import RESULTS_DIR, LEAKAGE_ANALYSIS_DIR

EXP_DIR     = os.path.join(RESULTS_DIR, "experiments")
LEAKAGE_DIR = LEAKAGE_ANALYSIS_DIR
CSV_DIR     = os.path.join(RESULTS_DIR, "csv")
os.makedirs(CSV_DIR, exist_ok=True)


def _safe_load(path):
    try:
        return json.load(open(path))
    except Exception:
        return None


def figure2():
    rows = []
    for path in sorted(glob(os.path.join(EXP_DIR, "figure2_*.json"))):
        d = _safe_load(path)
        if not d:
            continue
        label = d["label"]; n = d["n_monomers"]; t = d["t"]
        for k_str, rec in d.get("results", {}).items():
            k = int(k_str)
            rows.append(dict(
                system=label,
                n_monomers=n,
                t=t,
                k=k,
                num_blocks=rec.get("num_blocks"),
                projected_total_s=rec.get("projected_total"),
                probe_count=rec.get("probe_count"),
                probe_total_s=rec.get("probe_total"),
                wall_s=rec.get("wall"),
                error=rec.get("error", ""),
            ))
    _write(rows, "figure2_probe_curves.csv",
           ["system", "n_monomers", "t", "k", "num_blocks",
            "projected_total_s", "probe_count", "probe_total_s", "wall_s",
            "error"])


def figures3_4():
    p2 = _safe_load(os.path.join(EXP_DIR, "phase2_covering.json"))
    p3 = _safe_load(os.path.join(EXP_DIR, "phase3_full_hb.json"))

    rows = []
    if p2:
        for r in p2.get("results", []):
            full = r.get("full") or {}
            nt   = full.get("normaliz_time")
            bp   = r.get("best_projected")
            skip = r.get("full_run_skipped_reason")
            if nt is not None and nt > 0:
                value_s = nt
                source  = "full_run"
            elif bp is not None and bp not in (float("inf"),):
                value_s = bp
                source  = "probe_projected"
            else:
                value_s = None
                source  = "n/a"
            rows.append(dict(
                family=r["family"],
                size=r["size"],
                t=r["t"],
                n_monomers=r.get("n_monomers"),
                best_k=r.get("best_k"),
                value_s=value_s,
                source=source,
                full_normaliz_s=nt,
                full_wall_s=full.get("wall"),
                full_unique_vectors=full.get("unique_vectors"),
                best_projected_s=bp,
                skipped_reason=skip or "",
            ))
    _write(rows, "figure3_4_covering_bars.csv",
           ["family", "size", "t", "n_monomers", "best_k",
            "value_s", "source",
            "full_normaliz_s", "full_wall_s", "full_unique_vectors",
            "best_projected_s", "skipped_reason"])

    hb_rows = []
    if p3:
        for r in p3.get("results", []):
            hb_rows.append(dict(
                family=r["family"],
                size=r["size"],
                n_monomers=r.get("n_monomers"),
                normaliz_s=r.get("normaliz_time"),
                wall_s=r.get("wall"),
                vectors=r.get("vectors_found"),
                truncated=r.get("truncated"),
                hb_path=r.get("hb_path", ""),
            ))
    _write(hb_rows, "figure3_4_full_hb.csv",
           ["family", "size", "n_monomers",
            "normaliz_s", "wall_s", "vectors", "truncated", "hb_path"])


def leakage_vs_t():
    summary = _safe_load(os.path.join(LEAKAGE_DIR, "vary_t", "n7_incomplete", "summary.json"))
    if summary is None:
        print("[leakage_vs_t] missing summary.json -- skip"); return

    cutoff = summary.get("significance_cutoff_M")
    n_exp  = summary.get("n_expected_polymers")
    agg_rows = []
    for tag, agg in summary.get("by_set", {}).items():
        agg_rows.append(dict(
            set=tag,
            n_expected=n_exp,
            cutoff_M=cutoff,
            total_unexpected_conc_M=agg.get("total_unexpected_conc"),
            total_expected_deficit_M=agg.get("total_expected_deficit"),
            total_abs_deviation_M=agg.get("total_abs_deviation"),
            expected_recovered=agg.get("expected_recovered"),
            expected_missing=agg.get("expected_missing"),
            n_unexpected_above_cutoff=agg.get("n_unexpected_above"),
            n_unexpected_below_cutoff=agg.get("n_unexpected_below"),
        ))
    _write(agg_rows, "leakage_vs_t.csv",
           ["set", "n_expected", "cutoff_M",
            "total_unexpected_conc_M", "total_expected_deficit_M",
            "total_abs_deviation_M",
            "expected_recovered", "expected_missing",
            "n_unexpected_above_cutoff", "n_unexpected_below_cutoff"])

    long_rows = []
    base = os.path.join(LEAKAGE_DIR, "vary_t", "n7_incomplete")
    for f in sorted(glob(os.path.join(base, "polymer_compare_*.csv"))):
        tag = os.path.basename(f).replace("polymer_compare_", "").replace(".csv", "")
        with open(f) as fh:
            r = csv.DictReader(fh)
            for row in r:
                row["set"] = tag
                long_rows.append(row)
    if long_rows:
        cols = ["set"] + [c for c in long_rows[0] if c != "set"]
        _write(long_rows, "leakage_vs_t_per_polymer.csv", cols)


def leakage_vs_K():
    summary = _safe_load(os.path.join(LEAKAGE_DIR, "vary_removed_inputs", "n7_systems_compare",
                                      "summary.json"))
    if summary is None:
        print("[leakage_vs_K] missing summary.json -- skip"); return

    agg_rows = []
    for K_str, rec in summary.get("by_K", {}).items():
        agg = rec.get("aggregate", {})
        agg_rows.append(dict(
            removed_input_index=int(K_str),
            system_label=rec.get("system_label"),
            n_monomers=rec.get("n_monomers"),
            n_expected=rec.get("n_expected"),
            n_full_pstar=rec.get("n_full_pstar"),
            cutoff_M=agg.get("cutoff"),
            total_unexpected_conc_M=agg.get("total_unexpected_conc"),
            total_expected_deficit_M=agg.get("total_expected_deficit"),
            total_abs_deviation_M=agg.get("total_abs_deviation"),
            expected_recovered=agg.get("expected_recovered"),
            expected_missing=agg.get("expected_missing"),
            n_unexpected_above_cutoff=agg.get("n_unexpected_above"),
            n_unexpected_below_cutoff=agg.get("n_unexpected_below"),
        ))
    _write(agg_rows, "leakage_vs_K.csv",
           ["removed_input_index", "system_label", "n_monomers",
            "n_expected", "n_full_pstar", "cutoff_M",
            "total_unexpected_conc_M", "total_expected_deficit_M",
            "total_abs_deviation_M",
            "expected_recovered", "expected_missing",
            "n_unexpected_above_cutoff", "n_unexpected_below_cutoff"])

    base = os.path.join(LEAKAGE_DIR, "vary_removed_inputs", "n7_systems_compare")
    long_rows = []
    for K_dir in sorted(glob(os.path.join(base, "K*"))):
        K = os.path.basename(K_dir).lstrip("K")
        f = os.path.join(K_dir, "polymers_sorted.csv")
        if not os.path.exists(f):
            continue
        with open(f) as fh:
            r = csv.DictReader(fh)
            for row in r:
                row["removed_input_index"] = K
                long_rows.append(row)
    if long_rows:
        cols = ["removed_input_index"] + [c for c in long_rows[0]
                                       if c != "removed_input_index"]
        _write(long_rows, "leakage_vs_K_per_polymer.csv", cols)


def _write(rows, name, columns):
    path = os.path.join(CSV_DIR, name)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {path}  ({len(rows)} rows)")


def main():
    print(f"Writing CSVs to {CSV_DIR}\n")
    figure2()
    figures3_4()
    leakage_vs_t()
    leakage_vs_K()


if __name__ == "__main__":
    main()
