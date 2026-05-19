"""
make_plots.py — Render paper figures from whatever experiments have finished.

Reads:
  results/experiments/figure2_<system>.json   (Phase 1 probe curves)
  results/experiments/phase2_covering.json    (Phase 2 covering bars)
  results/experiments/phase3_full_hb.json     (Phase 3 Full-HB bars, if present)

Writes:
  results/figures/figure2_probe_curves.png    (Section 6.2 / paper Fig. 2)
  results/figures/figure3_linear_cascade.png  (Section 7.2 / paper Fig. 3)
  results/figures/figure4_binary_tree.png     (Section 7.2 / paper Fig. 4 left)
  results/figures/figure4_dna_cascade.png     (Section 7.2 / paper Fig. 4 right)

Each cell missing from JSON is silently skipped (so partial data still plots).
"""

import json
import os
import sys
from glob import glob


# Keep matplotlib cache inside results/ so scripts run cleanly on locked-down machines.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results", ".matplotlib"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # parent src/
from paths import RESULTS_DIR

EXP_DIR = os.path.join(RESULTS_DIR, "experiments")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Figure 2 — Probe-estimated total runtime vs k
# ---------------------------------------------------------------------------

def figure2_probe_curves():
    paths = sorted(glob(os.path.join(EXP_DIR, "figure2_*.json")))
    if not paths:
        print("[figure2] No figure2_*.json found — skipping.")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted_any = False

    for p in paths:
        data = json.load(open(p))
        label = data["label"]
        n     = data["n_monomers"]
        t     = data["t"]
        res   = data["results"]
        if not res:
            print(f"[figure2] {label}: empty results — skipping")
            continue

        ks, ys = [], []
        full_hb = None  # (k, y)
        for k_str, r in sorted(res.items(), key=lambda x: int(x[0])):
            k = int(k_str)
            if "error" in r:
                continue
            proj = r.get("projected_total")
            if proj is None:
                continue
            if k == n:
                full_hb = (k, proj)
            else:
                ks.append(k); ys.append(proj)

        if ks:
            line, = ax.plot(ks, ys, marker="o", linewidth=1.5, markersize=5,
                            label=f"{label}  (|M|={n})")
            plotted_any = True
            kmin = ks[np.argmin(ys)]
            ax.axvline(kmin, color=line.get_color(), linestyle=":",
                       alpha=0.4, linewidth=1)

        if full_hb is not None:
            ax.plot([full_hb[0]], [full_hb[1]], marker="*",
                    markersize=14,
                    label=f"{label}  Full HB  (k={full_hb[0]})")
            plotted_any = True

    if not plotted_any:
        plt.close(fig)
        print("[figure2] No usable data — no plot written.")
        return

    ax.set_xlabel("block size k")
    ax.set_ylabel("projected total runtime (s) — log scale")
    ax.set_yscale("log")
    ax.set_title(f"Probe-estimated Hilbert-basis runtime vs k  (t = 5)")
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    out = os.path.join(FIG_DIR, "figure2_probe_curves.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[figure2] wrote {out}")


# ---------------------------------------------------------------------------
# Figures 3 & 4 — grouped bar charts: covering bars per t + Full HB per size
# ---------------------------------------------------------------------------

T_VALUES_PLOT = [3, 4, 5, 6, 7, 8]
T_COLORS = plt.cm.viridis(np.linspace(0.15, 0.85, len(T_VALUES_PLOT)))
FULL_HB_COLOR = "#d62728"


def _load_phase_data():
    p2_path = os.path.join(EXP_DIR, "phase2_covering.json")
    p3_path = os.path.join(EXP_DIR, "phase3_full_hb.json")
    p2 = json.load(open(p2_path))["results"] if os.path.exists(p2_path) else []
    p3 = json.load(open(p3_path))["results"] if os.path.exists(p3_path) else []
    return p2, p3


def _cell_height(rec):
    """Return (height_seconds, is_projection_only) for a Phase-2 cell.

    Prefers the actual `full.normaliz_time` when the full enumeration ran;
    falls back to the probe-derived `best_projected` when the full run was
    skipped (projected > PHASE2_FULL_RUN_CAP_S). Returns (None, False) when
    no usable value exists.
    """
    full = rec.get("full") or {}
    nt = full.get("normaliz_time")
    if nt and nt > 0:
        return nt, False
    proj = rec.get("best_projected")
    if proj and proj > 0:
        return proj, True
    return None, False


def _bar_chart(family, sizes, title, out_name):
    p2, p3 = _load_phase_data()
    p2_fam = [r for r in p2 if r["family"] == family]
    p3_fam = {r["size"]: r for r in p3 if r["family"] == family}

    sizes_with_data = [s for s in sizes
                       if any(r["size"] == s for r in p2_fam) or s in p3_fam]
    if not sizes_with_data:
        print(f"[{out_name}] no data for family={family} — skipping")
        return

    fig, ax = plt.subplots(figsize=(max(7, 1.4 * len(sizes_with_data) + 3), 5))

    n_bars_per_group = len(T_VALUES_PLOT) + 1   # +1 for Full HB
    bar_w   = 0.85 / n_bars_per_group
    x_base  = np.arange(len(sizes_with_data))

    # Covering bars: real / projection-est / no-data
    used_projection_label = False
    used_missing_label    = False
    for i, t in enumerate(T_VALUES_PLOT):
        heights        = []
        is_projected_  = []
        is_missing_    = []           # cell ran but no covering was available
        for j, s in enumerate(sizes_with_data):
            rec = next((r for r in p2_fam if r["size"] == s and r["t"] == t), None)
            if rec is None:
                heights.append(np.nan); is_projected_.append(False); is_missing_.append(False); continue
            h, proj_only = _cell_height(rec)
            if h is not None:
                heights.append(h)
                is_projected_.append(proj_only)
                is_missing_.append(False)
            else:
                # Cell was attempted but no covering design existed for any k —
                # render as a small bar at the y-axis floor so the missing cell
                # is visible (rather than silently absent).
                heights.append(np.nan)
                is_projected_.append(False)
                is_missing_.append(True)
        offset = (i - n_bars_per_group / 2 + 0.5) * bar_w
        bars = ax.bar(x_base + offset, heights, bar_w,
                      label=f"t={t}", color=T_COLORS[i],
                      edgecolor="black", linewidth=0.3)
        # Hatch + small "est" annotation on projection-only bars.
        for bar_obj, proj_only, missing in zip(bars, is_projected_, is_missing_):
            if proj_only and not np.isnan(bar_obj.get_height()):
                bar_obj.set_hatch("///")
                bar_obj.set_alpha(0.55)
                ax.text(bar_obj.get_x() + bar_obj.get_width() / 2,
                        bar_obj.get_height() * 1.05, "est",
                        ha="center", va="bottom", fontsize=6,
                        color="dimgray", rotation=90)
                used_projection_label = True
            elif missing:
                # Mark missing slot with a hatched grey placeholder at the
                # bottom of the axis so the gap is unambiguous on inspection.
                ymin, _ = ax.get_ylim()
                xpos = bar_obj.get_x() + bar_obj.get_width() / 2
                ax.text(xpos, 1.05, "n/a",
                        ha="center", va="bottom", fontsize=6,
                        color="dimgray", rotation=90)
                used_missing_label = True

    # Full HB bars
    hb_heights, hb_truncated = [], []
    for s in sizes_with_data:
        rec = p3_fam.get(s)
        if rec is None:
            hb_heights.append(np.nan); hb_truncated.append(False)
        else:
            nt = rec.get("normaliz_time")
            hb_heights.append(nt if nt and nt > 0 else np.nan)
            hb_truncated.append(bool(rec.get("truncated", False)))
    offset = (n_bars_per_group / 2 - 0.5) * bar_w
    bars = ax.bar(x_base + offset, hb_heights, bar_w,
                  label="Full HB", color=FULL_HB_COLOR,
                  edgecolor="black", linewidth=0.3)
    # Mark truncated Full-HB bars with "hrs+"
    for j, (bar, truncated) in enumerate(zip(bars, hb_truncated)):
        if truncated and not np.isnan(bar.get_height()):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() * 1.08, "hrs+",
                    ha="center", va="bottom", fontsize=8, color=FULL_HB_COLOR)

    ax.set_xticks(x_base)
    ax.set_xticklabels([str(s) for s in sizes_with_data])
    ax.set_xlabel("system size")
    ax.set_ylabel("Normaliz time (s) — log scale")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.grid(True, which="both", axis="y", linestyle=":", alpha=0.4)

    notes = []
    if used_projection_label:
        notes.append("hatched / 'est': probe-projected time (full run skipped: projected > 30 min cap)")
    if used_missing_label:
        notes.append("'n/a': no LJCR covering design and DP fallback was off when Phase 2 ran")
    if notes:
        ax.text(0.02, 0.02, "\n".join(notes),
                transform=ax.transAxes, fontsize=8, color="dimgray",
                bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.8))

    fig.tight_layout()
    out = os.path.join(FIG_DIR, out_name)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[{out_name}] wrote {out}")


def figure3_linear_cascade():
    _bar_chart("cascade", [5, 6, 7, 8, 9],
               "Linear cascade — Pareto-optimal polymer enumeration time",
               "figure3_linear_cascade.png")


def figure4_binary_tree():
    _bar_chart("binary", [3, 4],
               "Binary tree — Pareto-optimal polymer enumeration time",
               "figure4_binary_tree.png")


def figure4_dna_cascade():
    _bar_chart("dna", [4, 5, 6, 7],
               "DNA cascade (2-2 modules) — Pareto-optimal polymer enumeration time",
               "figure4_dna_cascade.png")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    print(f"Writing figures to {FIG_DIR}\n")
    figure2_probe_curves()
    figure3_linear_cascade()
    figure4_binary_tree()
    figure4_dna_cascade()


if __name__ == "__main__":
    main()
