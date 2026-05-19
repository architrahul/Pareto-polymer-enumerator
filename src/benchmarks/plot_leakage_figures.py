"""
plot_leakage_figures.py — Reproduce paper Figures 5 and 6.

Compares each reduced candidate set's equilibrium concentrations against
the full Pareto-optimal set baseline, per polymer significant at
equilibrium. Plots relative concentration error (log scale) vs polymer
index (sorted by descending full-set concentration).

Expects the COFFEE outputs produced by leakage_coffee.py at
  results/leakage/coffee/n{N}_incomplete/{k25_t3, k25_t5, full_pstar}/

Writes:
  results/figures/figure5_leakage_t3.png        relative error using P̂_{25,3}
  results/figures/figure6_leakage_t5.png        relative error using P̂_{25,5}
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from paths import RESULTS_DIR, LEAKAGE_COFFEE_DIR

LEAKAGE_OUT = LEAKAGE_COFFEE_DIR
FIG_DIR     = os.path.join(RESULTS_DIR, "figures")

# Concentration threshold for "equilibrium-relevant" polymers (1 nM).
# Matches paper Section 7.3 footnote 7.
SIGNIFICANCE_THRESHOLD = 1e-9


def _load_polymers_and_concs(set_dir: str):
    """Return (polymer_vectors[N x M], concentrations[N])."""
    ocx_path = os.path.join(set_dir, "input.ocx")
    cof_path = os.path.join(set_dir, "coffee_output.txt")
    assert os.path.exists(ocx_path), ocx_path
    assert os.path.exists(cof_path), cof_path

    raw = np.loadtxt(ocx_path)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    # input.ocx columns: [index, 1, monomer counts..., energy]
    polymers = raw[:, 2:-1].astype(int)
    N        = polymers.shape[0]

    # COFFEE output: first line of N numbers = equilibrium concentrations
    concs = None
    with open(cof_path) as f:
        for line in f:
            parts = line.strip().split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) == N:
                concs = np.array(vals)
                break
    if concs is None:
        raise RuntimeError(f"could not parse concentrations from {cof_path} "
                           f"(expected a line of {N} floats)")
    return polymers, concs


def _per_polymer_relative_error(target_polys, target_concs,
                                reduced_polys, reduced_concs,
                                threshold=SIGNIFICANCE_THRESHOLD):
    """For each equilibrium-relevant polymer in the target (full P*) set,
    find its match in the reduced set and compute |Δc| / c_target."""
    mask = target_concs >= threshold
    sig_polys = target_polys[mask]
    sig_concs = target_concs[mask]
    order     = np.argsort(-sig_concs)
    sig_polys = sig_polys[order]
    sig_concs = sig_concs[order]

    # Build a lookup from polymer tuple → concentration in reduced set
    reduced_lookup = {tuple(p): c for p, c in zip(reduced_polys, reduced_concs)}

    rel_err = np.empty(len(sig_polys), dtype=float)
    present = np.empty(len(sig_polys), dtype=bool)
    for i, p in enumerate(sig_polys):
        c_red = reduced_lookup.get(tuple(p))
        if c_red is None:
            rel_err[i] = np.inf
            present[i] = False
        else:
            rel_err[i] = abs(c_red - sig_concs[i]) / sig_concs[i]
            present[i] = True
    return sig_polys, sig_concs, rel_err, present


def _plot_one(sig_concs, rel_err, present, title, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    n = len(sig_concs)
    x = np.arange(n)

    finite = present & np.isfinite(rel_err)
    # Plot present polymers as blue dots
    ax.scatter(x[finite], np.clip(rel_err[finite], 1e-16, None),
               s=14, color="#1f77b4", label="present in reduced set")

    # Plot absent polymers as red vertical bars at the top
    absent = ~present
    if absent.any():
        ax.axhline(1.0, color="lightgray", linestyle=":", linewidth=0.6)
        ymax_for_red = max(np.nanmax(rel_err[finite]) if finite.any() else 1.0, 1.0) * 10
        for xi in x[absent]:
            ax.plot([xi, xi], [SIGNIFICANCE_THRESHOLD, ymax_for_red],
                    color="#d62728", linewidth=0.8, alpha=0.8)
        # Single legend entry for absent polymers
        ax.plot([], [], color="#d62728", linewidth=2,
                label=f"absent in reduced set ({absent.sum()})")

    ax.set_yscale("log")
    ax.set_xlabel("polymer index (sorted by decreasing full-set concentration)")
    ax.set_ylabel("relative concentration error  |Δc| / c_target")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend(loc="best", fontsize=9)
    recovered = present.sum()
    ax.text(0.02, 0.95,
            f"{recovered}/{n} equilibrium-relevant polymers recovered "
            f"({recovered/n*100:.1f}%)",
            transform=ax.transAxes, fontsize=9,
            va="top", bbox=dict(boxstyle="round", fc="white", ec="0.8"))

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7.")
    args = parser.parse_args()

    base = os.path.join(LEAKAGE_OUT, f"n{args.cascade_n}_incomplete")
    print(f"Reading from {base}")
    os.makedirs(FIG_DIR, exist_ok=True)

    target_dir = os.path.join(base, "full_pstar")
    target_polys, target_concs = _load_polymers_and_concs(target_dir)
    print(f"  full P* : {len(target_polys)} polymers, "
          f"{(target_concs >= SIGNIFICANCE_THRESHOLD).sum()} above 1 nM")

    for set_label, fig_name, fig_caption in [
        ("k25_t3", "figure5_leakage_t3.png",
         f"Figure 5 — relative error using P̂_{{k=25, t=3}} on cascade-{args.cascade_n} leakage"),
        ("k25_t5", "figure6_leakage_t5.png",
         f"Figure 6 — relative error using P̂_{{k=25, t=5}} on cascade-{args.cascade_n} leakage"),
    ]:
        red_dir = os.path.join(base, set_label)
        if not os.path.exists(red_dir):
            print(f"  [{set_label}] missing — skip ({red_dir})")
            continue
        reduced_polys, reduced_concs = _load_polymers_and_concs(red_dir)
        sig_polys, sig_concs, rel_err, present = _per_polymer_relative_error(
            target_polys, target_concs, reduced_polys, reduced_concs
        )
        out_path = os.path.join(FIG_DIR, fig_name)
        _plot_one(sig_concs, rel_err, present, fig_caption, out_path)


if __name__ == "__main__":
    main()
