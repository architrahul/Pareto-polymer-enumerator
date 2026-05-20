#!/usr/bin/env python3
"""
Consolidated benchmark runner for the paper experiments.

This is the single entry point for the experiments described in the README.
It writes organized CSVs and figures under

    results/benchmarks/

Experiments:
  1. runtime-vs-k probe curves for cascade-7 and damien-10
  2. runtime-vs-t, best k by probing, plus Full-HB baselines
  3. equilibrium recovery: full P* vs t=3 and t=5 on cascade-7 incomplete
  4. leakage analysis: removed-input sweep and t sweep on cascade-7 incomplete

The script is resumable: most outputs are cached by the helper scripts or by
checking result CSVs. Delete the corresponding subdirectory under
results/benchmarks/ if you want a fully fresh run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results", ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results", ".cache"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["XDG_CACHE_HOME"], exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from paths import EXAMPLE_TBNS_DIR, RESULTS_DIR, HILBERT_BASIS_DIR, COFFEE_RESULTS_DIR
from hilbert_pipeline import (
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

BENCH_DIR = Path(RESULTS_DIR) / "benchmarks"
HB_DIR = Path(HILBERT_BASIS_DIR)
LOG_DIR = BENCH_DIR / "logs"
SCRIPT_DIR = Path(__file__).resolve().parent

T_VALUES = list(range(3, 9))
PHASE2_FULL_RUN_CAP_S = 30 * 60
EXP1_T = 5
CUTOFF_M = 1e-9
ENABLE_LOGS = False


class NullLog:
    def write(self, *_args, **_kwargs):
        return 0
    def flush(self):
        return None
    def __enter__(self):
        return self
    def __exit__(self, *_exc):
        return False


def open_benchmark_log(path: Path, header: str):
    if not ENABLE_LOGS:
        return NullLog()
    path.parent.mkdir(parents=True, exist_ok=True)
    f = path.open("w")
    f.write(header + "\n")
    f.flush()
    return f


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def ensure_dirs() -> None:
    for d in [BENCH_DIR, HB_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    if ENABLE_LOGS:
        LOG_DIR.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {path} ({len(rows)} rows)")


def append_jsonl(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(json.dumps(record, default=str) + "\n")


def monomer_path(family: str, size: int) -> str:
    if family == "linear_cascade":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{size}.txt")
    if family == "binary_tree":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_binary_tree_d{size}.txt")
    if family == "dna_cascade":
        return os.path.join(EXAMPLE_TBNS_DIR, f"monomers_dna_tbn_depth{size}.txt")
    raise ValueError(f"unknown family: {family}")


def load_system(path: str):
    monomers = load_monomers(path)
    domains = get_all_unique_domains(monomers)
    return monomers, domains


def k_values_runtime_vs_k(n: int, t: int, allow_after_25: bool) -> list[int]:
    vals = list(range(t + 1, min(25, n) + 1))
    if allow_after_25 and n > 25:
        vals.extend(range(30, n + 1, 5))
        # For cascade-7, n=57. The requested post-25 sweep is increments of 5,
        # but k=55 is too close to the full system and should be skipped.
        vals = [k for k in vals if k != 55]
        if n not in vals:
            vals.append(n)
    elif n not in vals and n <= 25:
        vals.append(n)
    return sorted(set(vals))


def run_subprocess(args: list[str]) -> None:
    args = list(args)
    if ENABLE_LOGS and args and args[0] in {"leakage_compute_all.py", "leakage_experiment_t.py", "leakage_analysis.py"}:
        args.append("--logs")
    print("$", " ".join(args))
    subprocess.run([sys.executable, "-u"] + args, cwd=SCRIPT_DIR, check=True)


def maybe_copy(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        print(f"copied {src} -> {dst}")


# ---------------------------------------------------------------------------
# Experiment 1 — runtime vs k
# ---------------------------------------------------------------------------

def experiment1_runtime_vs_k() -> None:
    out_dir = BENCH_DIR / "01_runtime_vs_k"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"01_runtime_vs_k_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    systems = [
        dict(system="linear_cascade_n7", monomer_file=os.path.join(EXAMPLE_TBNS_DIR, "monomers_cascade_n7.txt"), fallback_dp=True, after25=True),
        dict(system="damien_n10", monomer_file=os.path.join(EXAMPLE_TBNS_DIR, "monomers_damien_n10.txt"), fallback_dp=False, after25=False),
    ]

    rows: list[dict] = []
    start_input_listener()
    with open_benchmark_log(log_path, f"experiment1_runtime_vs_k started {datetime.now().isoformat()}") as log:
        for spec in systems:
            monomers, domains = load_system(spec["monomer_file"])
            n = len(monomers)
            for k in k_values_runtime_vs_k(n, EXP1_T, spec["after25"]):
                print(f"\n[exp1] {spec['system']} |M|={n} t={EXP1_T} k={k}")
                rec = dict(
                    experiment="runtime_vs_k",
                    system=spec["system"],
                    monomer_file=spec["monomer_file"],
                    n_monomers=n,
                    n_domains=len(domains),
                    t=EXP1_T,
                    k=k,
                    fallback_dp=spec["fallback_dp"],
                )
                t0 = time.time()
                try:
                    blocks = load_covering_blocks(n, k, EXP1_T, fallback_dp=spec["fallback_dp"])
                    projected, probe_times, num_blocks = probe_k(
                        k, EXP1_T, blocks, n, monomers, "monomer", domains, n,
                        best_projected=None, log=log,
                    )
                    rec.update(
                        num_blocks=num_blocks,
                        projected_total_s=projected,
                        probe_iterations=len(probe_times),
                        probe_normaliz_s=sum(probe_times),
                        probe_wall_s=time.time() - t0,
                        error="",
                    )
                except Exception as e:
                    rec.update(error=str(e), probe_wall_s=time.time() - t0)
                rows.append(rec)
                write_csv(out_dir / "runtime_vs_k.csv", rows, [
                    "experiment", "system", "monomer_file", "n_monomers", "n_domains",
                    "t", "k", "fallback_dp", "num_blocks", "projected_total_s",
                    "probe_iterations", "probe_normaliz_s", "probe_wall_s", "error",
                ])
                append_jsonl(out_dir / "runtime_vs_k.jsonl", rec)

    plot_runtime_vs_k(out_dir / "runtime_vs_k.csv", out_dir)


def plot_runtime_vs_k(csv_path: Path, out_dir: Path) -> None:
    by_system: dict[str, list[dict]] = {}
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row.get("projected_total_s"):
                by_system.setdefault(row["system"], []).append(row)

    for system, rows in by_system.items():
        rows = sorted(rows, key=lambda r: int(r["k"]))
        x = [int(r["k"]) for r in rows]
        y = [float(r["projected_total_s"]) for r in rows]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(x, y, marker="o", linewidth=1.5)
        ax.set_yscale("log")
        ax.set_xlabel("block size k")
        ax.set_ylabel("probe-estimated total runtime (s)")
        ax.set_title(f"Runtime vs k — {system}, t={EXP1_T}")
        ax.grid(True, which="both", linestyle=":", alpha=0.4)
        fig.tight_layout()
        out = out_dir / f"runtime_vs_k_{system}.png"
        fig.savefig(out, dpi=160)
        plt.close(fig)
        print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Experiment 2 — runtime vs t with best k + Full HB baselines
# ---------------------------------------------------------------------------

def experiment2_runtime_by_t() -> None:
    out_dir = BENCH_DIR / "02_runtime_by_t"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"02_runtime_by_t_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    targets = []
    targets += [dict(family="linear_cascade", size=s, fallback_dp=False) for s in range(5, 10)]
    targets += [dict(family="binary_tree", size=s, fallback_dp=True) for s in [3, 4]]
    targets += [dict(family="dna_cascade", size=s, fallback_dp=False) for s in range(4, 8)]

    covering_rows: list[dict] = []
    probe_rows: list[dict] = []
    full_rows: list[dict] = []

    start_input_listener()
    with open_benchmark_log(log_path, f"experiment2_runtime_by_t started {datetime.now().isoformat()}") as log:
        for target in targets:
            family, size = target["family"], target["size"]
            path = monomer_path(family, size)
            monomers, domains = load_system(path)
            n = len(monomers)

            for t in T_VALUES:
                print(f"\n[exp2] {family} size={size} |M|={n} t={t}")
                best = None
                best_blocks = None
                probe_wall_start = time.time()
                probe_normaliz_total = 0.0
                for k in range(t + 1, min(25, n) + 1):
                    rec = dict(family=family, size=size, n_monomers=n, t=t, k=k, fallback_dp=target["fallback_dp"])
                    t0 = time.time()
                    try:
                        blocks = load_covering_blocks(n, k, t, fallback_dp=target["fallback_dp"])
                        projected, probe_times, num_blocks = probe_k(
                            k, t, blocks, n, monomers, "monomer", domains, n,
                            best_projected=(best or {}).get("best_projected_s"), log=log,
                        )
                        probe_normaliz_total += sum(probe_times)
                        rec.update(
                            num_blocks=num_blocks,
                            projected_total_s=projected,
                            probe_iterations=len(probe_times),
                            probe_normaliz_s=sum(probe_times),
                            probe_wall_s=time.time() - t0,
                            pruned=projected is None,
                            error="",
                        )
                        if projected is not None and (best is None or projected < best["best_projected_s"]):
                            best = dict(best_k=k, best_projected_s=projected, num_blocks=num_blocks)
                            best_blocks = blocks
                    except Exception as e:
                        rec.update(error=str(e), probe_wall_s=time.time() - t0)
                    probe_rows.append(rec)
                    write_csv(out_dir / "probe_details.csv", probe_rows, [
                        "family", "size", "n_monomers", "t", "k", "fallback_dp", "num_blocks",
                        "projected_total_s", "probe_iterations", "probe_normaliz_s",
                        "probe_wall_s", "pruned", "error",
                    ])

                cell = dict(
                    family=family, size=size, monomer_file=path, n_monomers=n, n_domains=len(domains),
                    t=t, fallback_dp=target["fallback_dp"], best_k=None, best_projected_s=None,
                    best_num_blocks=None, total_probe_normaliz_s=probe_normaliz_total,
                    total_probe_wall_s=time.time() - probe_wall_start,
                    run_type="no_valid_k", actual_covering_wall_s=None,
                    actual_covering_normaliz_s=None, unique_vectors=None, hb_path="",
                    skipped_reason="no_valid_k", error="",
                )

                if best is not None:
                    cell.update(best_k=best["best_k"], best_projected_s=best["best_projected_s"], best_num_blocks=best["num_blocks"])
                    if best["best_projected_s"] <= PHASE2_FULL_RUN_CAP_S:
                        print(f"  best k={best['best_k']} projected={best['best_projected_s']:.1f}s <= 1800s; running full covering")
                        try:
                            result, _ = full_run_k(best["best_k"], best_blocks, n, monomers, "monomer", domains, n, log)
                            if result is not None:
                                hb_path = HB_DIR / f"exp2_{family}_size{size}_t{t}_k{best['best_k']}.txt"
                                save_polymer_vectors(result["vectors"], str(hb_path), n_monomers=n,
                                                     comment=f"exp2 covering {family} size={size} t={t} k={best['best_k']}")
                                cell.update(
                                    run_type="actual_covering",
                                    actual_covering_wall_s=result["total_wall_time"],
                                    actual_covering_normaliz_s=result["total_normaliz_time"],
                                    unique_vectors=result["unique_vectors"],
                                    hb_path=str(hb_path),
                                    skipped_reason="",
                                )
                        except Exception as e:
                            cell.update(run_type="actual_covering_failed", error=str(e))
                    else:
                        cell.update(run_type="estimated", skipped_reason="best_projected_over_1800s")
                        print(f"  best k={best['best_k']} projected={best['best_projected_s']:.1f}s > 1800s; using estimate")

                covering_rows.append(cell)
                write_csv(out_dir / "runtime_by_t_covering.csv", covering_rows, [
                    "family", "size", "monomer_file", "n_monomers", "n_domains", "t", "fallback_dp",
                    "best_k", "best_projected_s", "best_num_blocks", "total_probe_normaliz_s",
                    "total_probe_wall_s", "run_type", "actual_covering_wall_s",
                    "actual_covering_normaliz_s", "unique_vectors", "hb_path", "skipped_reason", "error",
                ])

            # Full-HB baseline: one Normaliz call on the complete monomer set.
            full_key = dict(family=family, size=size, monomer_file=path, n_monomers=n, n_domains=len(domains))
            print(f"\n[exp2/full] {family} size={size} |M|={n}")
            try:
                cleanup_normaliz_files()
                wall0 = time.time()
                elapsed, raw = run_normaliz_on_subset(monomers)
                wall = time.time() - wall0
                truncated = elapsed >= NORMALIZ_TIMEOUT_SECONDS - 1 or not raw
                hb_path = ""
                projected_vectors = []
                if raw and not truncated:
                    vecs = {tuple(v[:n]) for v in raw}
                    vecs.discard(tuple([0] * n))
                    hb = HB_DIR / f"exp2_full_hb_{family}_size{size}.txt"
                    save_polymer_vectors(vecs, str(hb), n_monomers=n,
                                         comment=f"exp2 Full HB {family} size={size}")
                    hb_path = str(hb)
                    projected_vectors = vecs
                full_rows.append(dict(**full_key, normaliz_s=elapsed, wall_s=wall,
                                      raw_vectors=len(raw), projected_vectors=len(projected_vectors),
                                      truncated=truncated, hb_path=hb_path, error=""))
            except Exception as e:
                full_rows.append(dict(**full_key, error=str(e)))
            write_csv(out_dir / "runtime_by_t_full_hb.csv", full_rows, [
                "family", "size", "monomer_file", "n_monomers", "n_domains", "normaliz_s",
                "wall_s", "raw_vectors", "projected_vectors", "truncated", "hb_path", "error",
            ])

    plot_runtime_by_t(out_dir / "runtime_by_t_covering.csv", out_dir / "runtime_by_t_full_hb.csv", out_dir)


def plot_runtime_by_t(cover_csv: Path, full_csv: Path, out_dir: Path) -> None:
    cover_rows = list(csv.DictReader(cover_csv.open())) if cover_csv.exists() else []
    full_rows = list(csv.DictReader(full_csv.open())) if full_csv.exists() else []

    for family in sorted({r["family"] for r in cover_rows}):
        rows = [r for r in cover_rows if r["family"] == family]
        sizes = sorted({int(r["size"]) for r in rows})
        fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(sizes)), 5))
        bar_w = 0.11
        group_x = list(range(len(sizes)))
        for i, t in enumerate(T_VALUES):
            vals = []
            hatches = []
            for size in sizes:
                match = next((r for r in rows if int(r["size"]) == size and int(r["t"]) == t), None)
                if not match:
                    vals.append(math.nan); hatches.append("")
                    continue
                if match.get("actual_covering_normaliz_s"):
                    vals.append(float(match["actual_covering_normaliz_s"])); hatches.append("")
                elif match.get("best_projected_s"):
                    vals.append(float(match["best_projected_s"])); hatches.append("//")
                else:
                    vals.append(math.nan); hatches.append("")
            xs = [x + (i - len(T_VALUES) / 2) * bar_w for x in group_x]
            bars = ax.bar(xs, vals, width=bar_w, label=f"t={t}", edgecolor="black", linewidth=0.25)
            for bar, hatch in zip(bars, hatches):
                bar.set_hatch(hatch)

        # Full HB as a black marker/bar at the right of each group.
        for x, size in zip(group_x, sizes):
            f = next((r for r in full_rows if r["family"] == family and int(r["size"]) == size), None)
            if f and f.get("normaliz_s"):
                ax.scatter([x + 0.42], [float(f["normaliz_s"])], marker="D", color="black", s=35,
                           label="Full HB" if x == group_x[0] else None)

        ax.set_yscale("log")
        ax.set_xticks(group_x)
        ax.set_xticklabels([str(s) for s in sizes])
        ax.set_xlabel("system size")
        ax.set_ylabel("runtime (s), log scale")
        ax.set_title(f"Runtime by t — {family}")
        ax.grid(True, which="both", axis="y", linestyle=":", alpha=0.4)
        ax.legend(fontsize=8, ncols=3)
        fig.tight_layout()
        out = out_dir / f"runtime_by_t_{family}.png"
        fig.savefig(out, dpi=160)
        plt.close(fig)
        print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Experiment 3 — equilibrium recovery
# ---------------------------------------------------------------------------

def experiment3_equilibrium_recovery() -> None:
    out_dir = BENCH_DIR / "03_equilibrium_recovery"
    figs_dir = out_dir / "figures"
    csv_dir = out_dir / "csv"
    for d in [figs_dir, csv_dir]:
        d.mkdir(parents=True, exist_ok=True)

    run_subprocess(["leakage_compute_all.py", "--n", "7"])
    run_subprocess(["leakage_coffee.py", "--n", "7"])
    run_subprocess(["plot_leakage_figures.py", "--n", "7"])
    run_subprocess(["dump_sorted_polymers.py", "--n", "7", "--top", "50"])

    maybe_copy(Path(RESULTS_DIR) / "figures" / "figure5_leakage_t3.png", figs_dir / "equilibrium_relative_error_t3.png")
    maybe_copy(Path(RESULTS_DIR) / "figures" / "figure6_leakage_t5.png", figs_dir / "equilibrium_relative_error_t5.png")

    # Keep large reusable files only in results/common/. This experiment folder
    # stores a compact index pointing to the shared Hilbert bases and COFFEE
    # outputs instead of copying hundreds of MB into a second location.
    hb_base = Path(HILBERT_BASIS_DIR)
    coffee_base = Path(COFFEE_RESULTS_DIR) / "n7_incomplete"
    rows = []
    for tag, hb_name in [
        ("full_pstar", "hilbert_full_p_star_n7_incomplete.txt"),
        ("k25_t3", "hilbert_k25_t3_monomer_n7_incomplete.txt"),
        ("k25_t5", "hilbert_k25_t5_monomer_n7_incomplete.txt"),
    ]:
        d = coffee_base / tag
        rows.append(dict(
            set=tag,
            shared_hilbert_basis=str(hb_base / hb_name),
            shared_coffee_dir=str(d),
            coffee_output=str(d / "coffee_output.txt"),
            sorted_csv=str(d / "coffee_output_sorted.csv"),
            polymers_sorted=str(d / "polymers_sorted.csv"),
        ))
    write_csv(csv_dir / "equilibrium_outputs.csv", rows, [
        "set", "shared_hilbert_basis", "shared_coffee_dir",
        "coffee_output", "sorted_csv", "polymers_sorted",
    ])

# ---------------------------------------------------------------------------
# Experiment 4 — leakage analysis
# ---------------------------------------------------------------------------

def experiment4_leakage() -> None:
    out_dir = BENCH_DIR / "04_leakage"
    removed_dir = out_dir / "removed_inputs"
    vary_t_dir = out_dir / "vary_t"
    csv_dir = out_dir / "csv"
    for d in [removed_dir, vary_t_dir, csv_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # 4.1 removed-input sweep K=1..7
    run_subprocess(["leakage_experiment_inputs.py", "--n", "7", "--K-values", "1", "2", "3", "4", "5", "6", "7"])
    src1 = Path(RESULTS_DIR) / "leakage" / "analysis" / "vary_removed_inputs" / "n7_systems_compare"
    maybe_copy(src1 / "summary.json", removed_dir / "summary.json")
    maybe_copy(src1 / "figure_leakage_vs_K.png", removed_dir / "leakage_vs_removed_inputs.png")
    for f in src1.glob("K*/polymers_sorted.csv"):
        maybe_copy(f, removed_dir / f.parent.name / "polymers_sorted.csv")
    for f in src1.glob("K*/polymer_compare.csv"):
        maybe_copy(f, removed_dir / f.parent.name / "polymer_compare.csv")

    # 4.2 t sweep for K=1 incomplete cascade
    run_subprocess(["leakage_experiment_t.py", "--n", "7"])
    src2 = Path(RESULTS_DIR) / "leakage" / "analysis" / "vary_t" / "n7_incomplete"
    maybe_copy(src2 / "summary.json", vary_t_dir / "summary.json")
    maybe_copy(src2 / "figure_leakage_vs_t.png", vary_t_dir / "leakage_vs_t.png")
    for f in src2.glob("polymer_compare*.csv"):
        maybe_copy(f, vary_t_dir / f.name)

    run_subprocess(["export_csv.py"])
    maybe_copy(Path(RESULTS_DIR) / "csv" / "leakage_vs_K.csv", csv_dir / "leakage_vs_removed_inputs.csv")
    maybe_copy(Path(RESULTS_DIR) / "csv" / "leakage_vs_t.csv", csv_dir / "leakage_vs_t.csv")
    maybe_copy(Path(RESULTS_DIR) / "csv" / "leakage_vs_K_per_polymer.csv", csv_dir / "leakage_vs_removed_inputs_per_polymer.csv")
    maybe_copy(Path(RESULTS_DIR) / "csv" / "leakage_vs_t_per_polymer.csv", csv_dir / "leakage_vs_t_per_polymer.csv")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--experiments", nargs="+", default=["1", "2", "3", "4"],
        choices=["1", "2", "3", "4", "runtime-vs-k", "runtime-by-t", "equilibrium", "leakage", "all"],
        help="Experiments to run. Default: 1 2 3 4.",
    )
    parser.add_argument("--logs", action="store_true", help="Write verbose benchmark/pipeline logs. Default: no verbose log files.")
    args = parser.parse_args()

    global ENABLE_LOGS
    ENABLE_LOGS = args.logs
    ensure_dirs()
    selected = set(args.experiments)
    if "all" in selected:
        selected = {"1", "2", "3", "4"}

    if "1" in selected or "runtime-vs-k" in selected:
        experiment1_runtime_vs_k()
    if "2" in selected or "runtime-by-t" in selected:
        experiment2_runtime_by_t()
    if "3" in selected or "equilibrium" in selected:
        experiment3_equilibrium_recovery()
    if "4" in selected or "leakage" in selected:
        experiment4_leakage()

    print(f"\nAll requested experiments complete. Organized outputs: {BENCH_DIR}")


if __name__ == "__main__":
    main()
