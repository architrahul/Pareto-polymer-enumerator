"""
generate_tbns.py

Generates four monomers.txt files for benchmarking:
  monomers_cascade_n<N>.txt   - linear cascade of depth N
  monomers_binary_tree_d<D>.txt - binary tree of depth D
  monomers_damien_n<N>.txt    - Damien 2-track sequential TBN of length N
  monomers_random_n<N>_m<M>.txt - random TBN with N domains, M monomers
"""

import random


# Cascade TBN

def generate_cascade(n: int) -> list[str]:
    """
    Linear cascade of depth n.
    Inputs: x1..x(n+1), each a pair of domains xi_1 xi_2.
    Module i computes y(i-1) + x(i+1) -> yi, for i in 1..n.
    y0 = x1 (first input signal).
    """
    monomers = []

    # Input signals
    for i in range(1, n + 2):
        monomers.append(f"x{i}_1 x{i}_2")

    def signal(name):
        return f"{name}_1", f"{name}_2"

    def module(a, b, c):
        a1, a2 = signal(a)
        b1, b2 = signal(b)
        c1, c2 = signal(c)
        return [
            f"{a1} {a2} {b1} {b2} {c1}",
            f"{a1}* {a2}* {b1}* {b2}*",
            f"{a2} {b1}",
            f"{b2} {c1} {c2}",
            f"{a2}* {b1}* {b2}* {c1}*",
            f"{c1}* {c2}*",
            f"{c1} {c2}",  # explicit output signal — paper Appendix D.1 "outputyi:"
        ]

    # Module 1: x1 + x2 -> y1
    monomers += module("x1", "x2", "y1")

    # Module i (i=2..n): y(i-1) + x(i+1) -> yi
    for i in range(2, n + 1):
        monomers += module(f"y{i-1}", f"x{i+1}", f"y{i}")

    return monomers

# DNA Strand Displacement TBNs

def generate_dna_tbn(depth: int) -> list[str]:
    monomers = []

    def signal(name):
        return f"{name}1", f"{name}2", f"{name}3"

    def module(a, b, c, d):
        a1, a2, a3 = signal(a)
        b1, b2, b3 = signal(b)
        c1, c2, c3 = signal(c)
        d1, d2, d3 = signal(d)

        return [
            # stage 1
            f"{a1} {a2} {a3} {b1} {b2} {b3} {c1} {d1}",
            f"{a1}* {a2}* {a3}* {b1}* {b2}* {b3}*",

            # stage 2
            f"{a2} {a3} {b2} {b3} {c1} {c2} {d1} {d2}",
            f"{a2}* {a3}* {b2}* {b3}* {c1}* {d1}*",

            # stage 3
            f"{a3} {b3} {c1} {c2} {c3} {d1} {d2} {d3}",
            f"{a3}* {b3}* {c1}* {c2}* {d1}* {d2}*",

            # sink layer
            f"{c1}* {c2}* {c3}* {d1}* {d2}* {d3}*",

            # explicit outputs
            f"{c1} {c2} {c3}",
            f"{d1} {d2} {d3}",
        ]

    # --- initial inputs (2 bundles) ---
    monomers.append(" ".join(signal("x1")))
    monomers.append(" ".join(signal("x2")))

    # --- first module ---
    monomers += module("x1", "x2", "y1", "z1")

    # --- cascade ---
    for i in range(2, depth + 1):
        monomers += module(
            f"y{i-1}", f"z{i-1}",
            f"y{i}", f"z{i}"
        )

    return monomers

# Binary Tree TBN

def generate_binary_tree(depth: int) -> list[str]:
    """
    Binary tree of given depth.
    Leaves are input signals. Each internal node combines two child signals
    into one parent signal using the same 6-monomer module as cascade.
    Total leaves = 2^depth, total internal nodes = 2^depth - 1.
    """
    monomers = []

    def signal(name):
        return f"{name}_1", f"{name}_2"

    def module(a, b, c):
        a1, a2 = signal(a)
        b1, b2 = signal(b)
        c1, c2 = signal(c)
        return [
            f"{a1} {a2} {b1} {b2} {c1}",
            f"{a1}* {a2}* {b1}* {b2}*",
            f"{a2} {b1}",
            f"{b2} {c1} {c2}",
            f"{a2}* {b1}* {b2}* {c1}*",
            f"{c1}* {c2}*",
            f"{c1} {c2}",  # explicit output signal — paper Appendix D.4 "outputyi:"
        ]

    n_leaves = 2 ** depth

    # Leaf input signals
    for i in range(1, n_leaves + 1):
        monomers.append(f"leaf{i}_1 leaf{i}_2")

    # Build tree bottom-up
    # Level 0 = leaves: leaf1..leaf(2^depth)
    # Level d = internal nodes combining level d-1 pairs
    current_level = [f"leaf{i}" for i in range(1, n_leaves + 1)]
    node_counter = 1

    while len(current_level) > 1:
        next_level = []
        for i in range(0, len(current_level), 2):
            a = current_level[i]
            b = current_level[i + 1]
            c = f"node{node_counter}"
            node_counter += 1
            monomers += module(a, b, c)
            next_level.append(c)
        current_level = next_level

    return monomers

# -------------------------
# Damien TBN
# -------------------------

def generate_damien(n: int) -> list[str]:
    """
    2-track sequential TBN of length n (from Damien's preprint).
    One substrate strand binds all selector domains.
    Two tracks (t0, t1) each propagate independently along the chain.
    Multiplicities from the CSV are dropped (pipeline doesn't use them).
    """
    monomers = []

    # Substrate: binds all selector complements
    monomers.append(" ".join(f"s_{i}*" for i in range(1, n + 1)))

    for i in range(1, n + 1):
        if i < n:
            monomers.append(f"t0_{i} s_{i} t0_{i+1}*")
            monomers.append(f"t1_{i} s_{i} t1_{i+1}*")
        else:
            monomers.append(f"t0_{i} s_{i}")
            monomers.append(f"t1_{i} s_{i}")

    return monomers


# -------------------------
# Random TBN
# -------------------------

def generate_random(n_domains: int, n_monomers: int,
                    min_sites: int = 2, max_sites: int = 5,
                    seed: int = 42) -> list[str]:
    """
    Random TBN with n_domains domain types and n_monomers monomers.
    Each monomer has between min_sites and max_sites binding sites,
    each chosen uniformly from {d, d*} for d in 1..n_domains.
    Duplicates are removed.
    """
    rng = random.Random(seed)
    domains = [f"d{i}" for i in range(1, n_domains + 1)]
    seen = set()
    monomers = []

    attempts = 0
    while len(monomers) < n_monomers and attempts < n_monomers * 100:
        attempts += 1
        k = rng.randint(min_sites, max_sites)
        sites = []
        for _ in range(k):
            d = rng.choice(domains)
            if rng.random() < 0.5:
                d = d + "*"
            sites.append(d)
        monomer = " ".join(sites)
        if monomer not in seen:
            seen.add(monomer)
            monomers.append(monomer)

    if len(monomers) < n_monomers:
        print(f"Warning: only generated {len(monomers)}/{n_monomers} unique monomers")

    return monomers


# -------------------------
# Writer
# -------------------------

def write_monomers(filename: str, monomers: list[str], description: str = ""):
    with open(filename, "w") as f:
        if description:
            f.write(f"# {description}\n")
        for m in monomers:
            f.write(m + "\n")
    print(f"Written {len(monomers)} monomers to {filename}")


# -------------------------
# Main
# -------------------------

if __name__ == "__main__":
    # Cascade: depths 2, 5..10. Each has 8n+1 monomers.
    for n in [2, 5, 6, 7, 8, 9, 10]:
        m = generate_cascade(n)
        write_monomers(f"monomers_cascade_n{n}.txt", m, f"Linear cascade depth {n}")

    # Cascade with the first input monomer removed — paper Section 7.3 leakage
    # analysis. Paper reports |M|=56 for the "m=8 leakage" case, but m=8
    # incomplete has 64 monomers; m=7 incomplete has 56. We generate both so
    # we can match timings against the paper to identify the typo.
    for n in [7, 8]:
        m = generate_cascade(n)
        incomplete = [line for line in m if line != "x1_1 x1_2"]
        write_monomers(
            f"monomers_cascade_n{n}_incomplete.txt", incomplete,
            f"Linear cascade depth {n} with first input (x1) removed — leakage analysis",
        )

    # DNA TBN: depths
    for d in [3, 4, 5, 6, 7, 8]:
        m = generate_dna_tbn(d)
        write_monomers(f"monomers_dna_tbn_depth{d}.txt", m, f"DNA TBN cascade depth {d}")

    # Binary tree: depths 2, 3, 4
    for d in [2, 3, 4]:
        m = generate_binary_tree(d)
        write_monomers(f"monomers_binary_tree_d{d}.txt", m, f"Binary tree depth {d}")

    # Damien: lengths 5, 10
    for n in [5, 10]:
        m = generate_damien(n)
        write_monomers(f"monomers_damien_n{n}.txt", m, f"Damien 2-track TBN length {n}")

