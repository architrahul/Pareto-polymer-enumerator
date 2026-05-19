#!/usr/bin/env python3
"""
coffee_pipeline.py

Generate COFFEE inputs from monomer/polymer files and run coffee-cli.

every monomer concentration is 1 µM and binding energies are -20 per bond

Example:

python src/coffee_pipeline.py \
  --monomers example-tbns/monomers_cascade_n7.txt \
  --polymers results/common/hilbert_basis/phase3_full_hb_cascade_size7.txt \
  --out-dir results/linear_cascade_size7/coffee_input/normal \
  --coffee-cli coffee/crates/coffee-cli/target/release/coffee-cli \
  --label normal
"""

import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coffee_parser import (
    assign_domain_energies,
    generate_coffee_inputs,
    modify_output,
    parse_monomers,
)


def run_coffee(coffee_cli: str, ocx_path: str, con_path: str, out_path: str):
    assert os.path.isfile(coffee_cli), f"coffee-cli not found at {coffee_cli}"

    t0 = time.time()
    rc = subprocess.run(
        [coffee_cli, ocx_path, con_path, "-o", out_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )

    print(
        f"coffee-cli rc={rc.returncode} in {time.time() - t0:.1f}s "
        f"→ {out_path}"
    )

    if rc.returncode != 0:
        raise RuntimeError("coffee-cli failed")


def main():
    parser = argparse.ArgumentParser(
        description="Generate COFFEE inputs from Hilbert-basis output and run COFFEE."
    )

    parser.add_argument("--monomers", required=True, help="Path to monomer file.")
    parser.add_argument("--polymers", required=True, help="Path to polymer vector file.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--coffee-cli", required=True, help="Path to coffee-cli binary.")
    parser.add_argument("--label", default="coffee", help="Log label.")

    args = parser.parse_args()


    monomers = parse_monomers(args.monomers)
    print(f"Parsed monomers: {len(monomers)}")

    with open(args.polymers) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            vec = list(map(int, line.split()))
            print(f"First polymer vector length: {len(vec)}")
            print(f"First polymer line number: {line_no}")
            break

    domain_energy = assign_domain_energies(monomers, seed=42)
    
    generate_coffee_inputs(
        monomers=monomers,
        polymer_file=args.polymers,
        out_dir=args.out_dir,
        domain_energy=domain_energy,
        zero_indices=None,
        label=args.label,
    )

    ocx_path = os.path.join(args.out_dir, "input.ocx")
    con_path = os.path.join(args.out_dir, "input.con")
    coffee_out = os.path.join(args.out_dir, "coffee_output.txt")

    run_coffee(
        coffee_cli=args.coffee_cli,
        ocx_path=ocx_path,
        con_path=con_path,
        out_path=coffee_out,
    )
    csv_out = modify_output(coffee_out, args.polymers)
    print(f"Wrote sorted concentration/vector CSV → {csv_out}")


if __name__ == "__main__":
    main()
