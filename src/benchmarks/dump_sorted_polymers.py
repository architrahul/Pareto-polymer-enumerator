"""
dump_sorted_polymers.py — Dump (polymer, equilibrium concentration) rows
sorted by decreasing concentration, for each of the three leakage polymer
sets (k=25 t=3, k=25 t=5, full P*) produced by leakage_coffee.py.

Outputs (one per set, under results/leakage/coffee/n{N}_incomplete/<set>/):
  polymers_sorted.csv        Full sorted table: rank,concentration,energy,m1,m2,...
  polymers_sorted_top50.txt  Human-readable top-50 with monomer-name expansion

Each CSV row:
  rank, concentration_M, energy_kT, count_of_m1, count_of_m2, ...
column headers include the monomer index. A separate `monomer_names.txt`
in the parent dir maps index → monomer string (so columns are decoded).
"""

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coffee_parser import parse_monomers
from paths import EXAMPLE_TBNS_DIR, LEAKAGE_COFFEE_DIR

LEAKAGE_OUT = LEAKAGE_COFFEE_DIR
SETS = ["k25_t3", "k25_t5", "full_pstar"]


def _load_polymers_and_concs(set_dir: str):
    """Read input.ocx + coffee_output.txt → (polymers[N,M], energies[N], concs[N])."""
    ocx_path = os.path.join(set_dir, "input.ocx")
    cof_path = os.path.join(set_dir, "coffee_output.txt")
    raw = np.loadtxt(ocx_path)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    polymers = raw[:, 2:-1].astype(int)
    energies = raw[:,   -1].astype(float)
    N = polymers.shape[0]

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
        raise RuntimeError(f"could not parse concentrations from {cof_path}")
    return polymers, energies, concs


def _monomer_string(monomer_list):
    """Render a monomer (list of domain strings) as space-separated."""
    return " ".join(monomer_list)


def dump_set(set_dir: str, monomer_names, top_n: int = 50):
    print(f"\n  [{os.path.basename(set_dir)}]")
    polymers, energies, concs = _load_polymers_and_concs(set_dir)
    M = polymers.shape[1]
    assert M == len(monomer_names), \
        f"polymer width {M} != monomer count {len(monomer_names)}"

    order = np.argsort(-concs)
    polymers = polymers[order]
    energies = energies[order]
    concs    = concs[order]

    csv_path = os.path.join(set_dir, "polymers_sorted.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["rank", "concentration_M", "energy_kT"] + [f"m{i}" for i in range(M)]
        w.writerow(header)
        for rank, (p, e, c) in enumerate(zip(polymers, energies, concs)):
            w.writerow([rank, f"{c:.6e}", f"{e:.6e}"] + list(p))
    print(f"    wrote {csv_path}  ({len(polymers)} rows)")

    top_path = os.path.join(set_dir, f"polymers_sorted_top{top_n}.txt")
    with open(top_path, "w") as f:
        f.write(f"# Top {top_n} polymers by equilibrium concentration\n")
        f.write(f"# total polymers in set: {len(polymers)}\n")
        f.write(f"# columns: rank | concentration | energy_kT | polymer (n x monomer_string)\n\n")
        for rank, (p, e, c) in enumerate(zip(polymers, energies, concs)):
            if rank >= top_n:
                break
            parts = []
            for i, count in enumerate(p):
                if count > 0:
                    parts.append(f"{int(count)}x ({monomer_names[i]})")
            polymer_str = " + ".join(parts) if parts else "(empty)"
            f.write(f"{rank:>4}  c={c:.3e} M  G={e:.2e} kT  :: {polymer_str}\n")
    print(f"    wrote {top_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7.")
    parser.add_argument("--top", type=int, default=50,
                        help="Top-N polymers for the human-readable txt. Default 50.")
    args = parser.parse_args()

    label    = f"n{args.cascade_n}_incomplete"
    base_dir = os.path.join(LEAKAGE_OUT, label)
    assert os.path.isdir(base_dir), f"missing {base_dir}"

    monomer_file = os.path.join(EXAMPLE_TBNS_DIR,
                                f"monomers_cascade_n{args.cascade_n}_incomplete.txt")
    monomers      = parse_monomers(monomer_file)
    monomer_names = [_monomer_string(m) for m in monomers]
    print(f"Cascade {label}: |M|={len(monomer_names)} monomers")

    name_path = os.path.join(base_dir, "monomer_names.txt")
    with open(name_path, "w") as f:
        f.write(f"# Monomer index -> string  (cascade {label})\n")
        for i, s in enumerate(monomer_names):
            f.write(f"{i}\t{s}\n")
    print(f"wrote {name_path}")

    for s in SETS:
        set_dir = os.path.join(base_dir, s)
        if not os.path.isdir(set_dir):
            print(f"\n  [{s}] missing - skip")
            continue
        dump_set(set_dir, monomer_names, top_n=args.top)


if __name__ == "__main__":
    main()
