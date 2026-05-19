import random
import re
import os
import csv


# -------------------------
# Parse monomer file
# -------------------------
def parse_monomers(path):
    monomers = []
    domain_pattern = re.compile(r"[a-zA-Z]+\d+_\d+\*?")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if ":" in line:
                line = line.split(":", 1)[1].strip()

            parts = line.split(",")
            domains = domain_pattern.findall(parts[0])
            monomers.append(domains)

    return monomers


# -------------------------
# Identify input signal monomers
# -------------------------
def find_input_signal_indices(monomers):
    """
    Input signal monomers in the cascade-style TBNs are the single-stranded
    input pairs `xN_1 xN_2` (exactly two domains, same xN prefix). Anything
    else — module monomers, intermediates like `x1_2 x2_1`, or anything with
    starred / y / z / leaf / node domains — is excluded.
    """
    pair_pattern = re.compile(r"^x(\d+)_(\d+)$")
    input_indices = []

    for i, m in enumerate(monomers):
        if len(m) != 2:
            continue
        matches = [pair_pattern.match(d) for d in m]
        if not all(matches):
            continue
        # Same xN prefix on both, and the two suffixes form the pair {1, 2}
        x_ids   = {mt.group(1) for mt in matches}
        suffix  = {mt.group(2) for mt in matches}
        if len(x_ids) == 1 and suffix == {"1", "2"}:
            input_indices.append(i)

    return input_indices


# --------------------
# Assign bond energies
# --------------------
def assign_domain_energies(monomers, seed=42):
    random.seed(seed)

    domain_types = set()
    for m in monomers:
        for d in m:
            base = d.replace("*", "")
            domain_types.add(base)

    energies = {}
    for d in sorted(domain_types):
        energies[d] = -20.0

    return energies


# -------------------------
# Count bonds in polymer
# -------------------------
def compute_polymer_energy(polymer, monomers, domain_energy):
    domain_count = {}

    for i, count in enumerate(polymer):
        if count == 0:
            continue
        for _ in range(count):
            for d in monomers[i]:
                domain_count[d] = domain_count.get(d, 0) + 1

    energy = 0.0
    for d in list(domain_count.keys()):
        if d.endswith("*"):
            base = d[:-1]
            if base in domain_count:
                bonds = min(domain_count[d], domain_count[base])
                energy += bonds * domain_energy[base]

    return energy


# -------------------------
# Read polymer file
# -------------------------
def read_polymers(path):
    polymers = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vec = list(map(int, line.split()))
            polymers.append(vec)

    return polymers


# -------------------------
# Convert raw COFFEE output into an analysis-friendly CSV
# -------------------------
def modify_output(coffee_output_path, polymer_file, csv_output_path=None):
    """Write concentrations paired with polymer vectors, sorted descending.

    COFFEE's raw output is kept untouched because the plotting scripts parse
    it directly. This companion CSV is easier to inspect manually:

        concentration_M,polymer_vector
        2.630000e-09,0 1 1 1 ...

    Returns the CSV path.
    """
    polymers = read_polymers(polymer_file)
    n_polymers = len(polymers)

    concs = None
    with open(coffee_output_path) as f:
        for line in f:
            parts = line.strip().split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) == n_polymers:
                concs = vals
                break

    if concs is None:
        raise RuntimeError(
            f"could not parse {n_polymers} concentrations from {coffee_output_path}"
        )

    if csv_output_path is None:
        csv_output_path = os.path.join(
            os.path.dirname(coffee_output_path),
            "coffee_output_sorted.csv",
        )

    rows = sorted(zip(concs, polymers), key=lambda item: item[0], reverse=True)
    with open(csv_output_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["concentration_M", "polymer_vector"])
        for conc, vec in rows:
            w.writerow([f"{conc:.6e}", " ".join(str(x) for x in vec)])

    return csv_output_path


# -------------------------
# Write input.ocx
# -------------------------
def write_ocx(polymers, energies, output_path):
    with open(output_path, "w") as f:
        for i, (vec, energy) in enumerate(zip(polymers, energies), start=1):
            row = [str(i), "1"] + [str(x) for x in vec] + [f"{energy:.6e}"]
            f.write("\t".join(row) + "\n")


# -------------------------
# Write input.con
# normal: all monomers at 1 µM
# leak:   input signal monomers at 0, everything else at 1 µM
# -------------------------
def write_con(n_monomers, output_path, zero_indices=None):
    """
    Writes a .con file.
    - zero_indices: set/list of monomer indices (0-based) whose concentration
                    should be 0 (used for the leak experiment).
    All other monomers are set to 1 µM.
    """
    zero_set = set(zero_indices) if zero_indices else set()

    with open(output_path, "w") as f:
        for i in range(n_monomers):
            if i in zero_set:
                f.write("0.000000e+00\n")
            else:
                f.write("1.000000e-06\n")


# -------------------------
# Optional: save energies for reference
# -------------------------
def write_domain_energies(domain_energy, output_path):
    with open(output_path, "w") as f:
        for d in sorted(domain_energy):
            f.write(f"{d}\t{domain_energy[d]:.6f}\n")


# -------------------------
# Main pipeline
# -------------------------
def generate_coffee_inputs(monomers, polymer_file, out_dir, domain_energy,
                           zero_indices=None, label=""):
    os.makedirs(out_dir, exist_ok=True)

    polymers = read_polymers(polymer_file)

    energies = [
        compute_polymer_energy(p, monomers, domain_energy)
        for p in polymers
    ]

    write_ocx(polymers, energies, os.path.join(out_dir, "input.ocx"))
    write_con(len(monomers), os.path.join(out_dir, "input.con"),
              zero_indices=zero_indices)
    write_domain_energies(domain_energy, os.path.join(out_dir, "domain_energies.txt"))

    n_zero = len(zero_indices) if zero_indices else 0
    print(f"[{label}] Generated COFFEE inputs in: {out_dir}")
    print(f"  Polymers : {len(polymers)}")
    print(f"  Monomers : {len(monomers)}  ({n_zero} set to 0 conc, "
          f"{len(monomers) - n_zero} at 1 µM)")


# -------------------------
# Usage
# -------------------------
if __name__ == "__main__":

    from paths import EXAMPLE_TBNS_DIR, RESULTS_DIR

    monomer_file = os.path.join(EXAMPLE_TBNS_DIR, "monomers_cascade_n10.txt")
    polymer_file = os.path.join(
        RESULTS_DIR, "linear_cascade_m10",
        "hilbert_81_k25_t5_covering.txt"
    )
    base_out_dir = os.path.join(
        RESULTS_DIR, "linear_cascade_m10", "coffee_input"
    )

    # Parse monomers and assign energies once
    monomers = parse_monomers(monomer_file)
    domain_energy = assign_domain_energies(monomers, seed=42)

    # Find which monomers are input signals (x1..x11 strands)
    input_signal_indices = find_input_signal_indices(monomers)
    print(f"Detected {len(input_signal_indices)} input signal monomers "
          f"(indices: {input_signal_indices})")

    # Normal run: all monomers at 1 µM
    generate_coffee_inputs(
        monomers=monomers,
        polymer_file=polymer_file,
        out_dir=os.path.join(base_out_dir, "normal"),
        domain_energy=domain_energy,
        zero_indices=None,
        label="normal",
    )

    # Leak run: input signal monomers at 0, everything else at 1 µM
    generate_coffee_inputs(
        monomers=monomers,
        polymer_file=polymer_file,
        out_dir=os.path.join(base_out_dir, "leak"),
        domain_energy=domain_energy,
        zero_indices=input_signal_indices,
        label="leak",
    )
