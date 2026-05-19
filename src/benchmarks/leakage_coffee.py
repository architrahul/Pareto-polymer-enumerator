"""
leakage_coffee.py — Build COFFEE inputs and run COFFEE on the three polymer
sets produced by leakage_compute_all.py:
  1. P̂_{k=25, t=3}  →  Figure 5 reduced candidate
  2. P̂_{k=25, t=5}  →  Figure 6 reduced candidate
  3. Full P*          →  baseline for both figures

Uses the *leakage* concentration scenario from paper Section 7.3: the first
input monomer is removed from the system, so it is absent from the incomplete
monomer file. Every remaining monomer, including surviving input monomers, is
set to 1 µM. Each polymer's free energy is computed by counting bonds at
-20 kT per bond, matching the paper's setup.

Outputs land in
  results/common/coffee/n{N}_incomplete/{set_label}/
    ├── input.ocx
    ├── input.con
    ├── domain_energies.txt
    └── coffee_output.txt   (after the COFFEE CLI runs)
"""

import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from paths import EXAMPLE_TBNS_DIR, LEAKAGE_HB_DIR, LEAKAGE_COFFEE_DIR, REPO_ROOT
from coffee_parser import (
    assign_domain_energies,
    generate_coffee_inputs,
    modify_output,
    parse_monomers,
)

COFFEE_CLI = os.path.join(REPO_ROOT, "coffee", "crates", "coffee-cli",
                          "target", "release", "coffee-cli")
LEAKAGE_OUT = LEAKAGE_COFFEE_DIR


def run_coffee(ocx_path: str, con_path: str, out_path: str):
    """Invoke coffee-cli with -o flag so out_path gets just the results line."""
    assert os.path.isfile(COFFEE_CLI), f"coffee-cli not built at {COFFEE_CLI}"
    t0 = time.time()
    rc = subprocess.run(
        [COFFEE_CLI, ocx_path, con_path, "-o", out_path],
        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    print(f"    coffee-cli rc={rc.returncode} in {time.time()-t0:.1f}s "
          f"→ {os.path.basename(out_path)}")


def setup_and_run(cascade_n: int):
    label = f"n{cascade_n}_incomplete"
    monomer_file = os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{cascade_n}_incomplete.txt")
    assert os.path.exists(monomer_file), f"missing {monomer_file}"

    polymer_files = {
        "k25_t3":    os.path.join(LEAKAGE_HB_DIR, f"hilbert_k25_t3_monomer_{label}.txt"),
        "k25_t5":    os.path.join(LEAKAGE_HB_DIR, f"hilbert_k25_t5_monomer_{label}.txt"),
        "full_pstar":os.path.join(LEAKAGE_HB_DIR, f"hilbert_full_p_star_{label}.txt"),
    }

    monomers      = parse_monomers(monomer_file)
    domain_energy = assign_domain_energies(monomers, seed=42)
    print(f"\nSystem: cascade {label}  |M|={len(monomers)}  "
          f"(first input absent; all remaining monomers at 1 µM)")

    base_dir = os.path.join(LEAKAGE_OUT, label)
    os.makedirs(base_dir, exist_ok=True)

    for set_label, p_path in polymer_files.items():
        if not os.path.exists(p_path):
            print(f"  [{set_label}] MISSING polymer file: {p_path} — SKIP")
            continue

        out_dir = os.path.join(base_dir, set_label)
        ocx     = os.path.join(out_dir, "input.ocx")
        con     = os.path.join(out_dir, "input.con")
        cof_out = os.path.join(out_dir, "coffee_output.txt")

        if os.path.exists(cof_out):
            csv_out = modify_output(cof_out, p_path)
            print(f"  [{set_label}] coffee_output already exists — "
                  f"wrote {os.path.basename(csv_out)}")
            continue

        print(f"\n  [{set_label}] building inputs from {os.path.basename(p_path)} ...")
        generate_coffee_inputs(
            monomers=monomers,
            polymer_file=p_path,
            out_dir=out_dir,
            domain_energy=domain_energy,
            zero_indices=None,
            label=set_label,
        )
        print(f"  [{set_label}] running coffee-cli ...")
        run_coffee(ocx, con, cof_out)
        csv_out = modify_output(cof_out, p_path)
        print(f"  [{set_label}] wrote {os.path.basename(csv_out)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=7, dest="cascade_n",
                        help="Cascade depth. Default 7 (paper's |M|=56 case).")
    args = parser.parse_args()
    setup_and_run(args.cascade_n)


if __name__ == "__main__":
    main()
