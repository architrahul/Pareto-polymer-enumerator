"""
expected_polymers.py — Build the "expected" polymer set + concentrations
for a linear-cascade system under leakage (some early inputs removed).

The leakage scenario: cascade with n modules and n+1 inputs (x1..x_{n+1}),
where the first `removed_count` inputs are dropped (concentration = 0).
With no first input(s), no cascade module fires its output, so each module
sits in its "leakage" equilibrium — the seven module monomers self-assemble
into three internal polymers. Plus the surviving inputs are each present as
free singletons.

Per-module groupings (m1..m7 per the paper / user's convention):
    m1 = a1* a2* b1* b2*            (the "complement" strand)
    m2 = a1  a2  b1  b2  c1         (the "central" strand)
    m3 = a2* b1* b2* c1*            (intermediate)
    m4 = a2  b1                     (intermediate)
    m5 = b2  c1  c2                 (intermediate)
    m6 = c1* c2*                    (output complement)
    m7 = c1  c2                     (output)

Expected polymers under leakage (per module):
    {m1, m2}           — central + its complement (4 a/b bonds, leaves c1)
    {m3, m4, m5}       — three intermediate strands (4 bonds, leaves c2 on m5)
    {m6, m7}           — output + its complement  (2 c bonds)

In `tbn_builder.py`'s output ordering, the monomers within a module are laid
out in this order (per `generate_cascade.module`):

    offset 0: m2  (central)
    offset 1: m1  (complement)
    offset 2: m4  (intermediate a-b)
    offset 3: m5  (intermediate b-c)
    offset 4: m3  (intermediate ax*-by*-c1*)
    offset 5: m6  (output complement)
    offset 6: m7  (output)

So within a module the three expected polymers occupy index offsets
{0, 1}, {2, 3, 4}, {5, 6}.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coffee_parser import parse_monomers
from paths import EXAMPLE_TBNS_DIR


# Within-module index offsets for each expected polymer (in tbn_builder layout):
MODULE_POLYMER_OFFSETS = [
    [0, 1],        # {m2, m1}     — central + complement
    [2, 3, 4],     # {m4, m5, m3} — intermediate triple
    [5, 6],        # {m6, m7}     — output complement + output
]
MONOMERS_PER_MODULE = 7

DEFAULT_INITIAL_CONC = 1e-6  # 1 µM — matches the existing leakage pipeline


@dataclass(frozen=True)
class ExpectedPolymer:
    label: str
    vector: tuple              # length = total monomer count, 1's at member indices
    concentration: float


def build_expected_polymers(
    cascade_n: int,
    removed_inputs: int = 1,
    initial_conc: float = DEFAULT_INITIAL_CONC,
) -> list[ExpectedPolymer]:
    """Return the list of expected polymers under leakage for cascade-n with the
    first `removed_inputs` input monomers (x1, x2, …, x_removed_inputs) removed.
    """
    if removed_inputs < 1 or removed_inputs > cascade_n + 1:
        raise ValueError(f"removed_inputs must be in [1, {cascade_n + 1}], got {removed_inputs}")

    fname = _monomer_filename(cascade_n, removed_inputs)
    path  = os.path.join(EXAMPLE_TBNS_DIR, fname)
    monomers = parse_monomers(path)
    n_total  = len(monomers)

    n_inputs_present = (cascade_n + 1) - removed_inputs
    expected_polys: list[ExpectedPolymer] = []

    # 1) Free-input singletons
    for i in range(n_inputs_present):
        vec = [0] * n_total
        vec[i] = 1
        expected_polys.append(ExpectedPolymer(
            label=f"input_x{removed_inputs + 1 + i}",
            vector=tuple(vec),
            concentration=initial_conc,
        ))

    # 2) Per-module internal polymers
    module_base = n_inputs_present  # first module monomer
    for mod_i in range(1, cascade_n + 1):
        module_start = module_base + (mod_i - 1) * MONOMERS_PER_MODULE
        for poly_idx, offsets in enumerate(MODULE_POLYMER_OFFSETS):
            vec = [0] * n_total
            for off in offsets:
                vec[module_start + off] = 1
            kind = ["central_complement", "intermediate", "output_pair"][poly_idx]
            expected_polys.append(ExpectedPolymer(
                label=f"module{mod_i}_{kind}",
                vector=tuple(vec),
                concentration=initial_conc,
            ))

    # Sanity check
    coverage = [0] * n_total
    for ep in expected_polys:
        for i, c in enumerate(ep.vector):
            coverage[i] += c
    missing = [i for i, c in enumerate(coverage) if c == 0]
    if missing:
        raise RuntimeError(
            f"expected polymer construction left {len(missing)} monomer(s) "
            f"uncovered (indices {missing[:5]}...) - module-offset map likely wrong"
        )

    return expected_polys


def _monomer_filename(cascade_n: int, removed_inputs: int) -> str:
    """Filename for cascade-n with first K inputs removed.

    K=1 keeps the existing  monomers_cascade_n{N}_incomplete.txt name (already
    in example-tbns/). K>=2 uses monomers_cascade_n{N}_incomplete{K}.txt.
    """
    if removed_inputs == 1:
        return f"monomers_cascade_n{cascade_n}_incomplete.txt"
    return f"monomers_cascade_n{cascade_n}_incomplete{removed_inputs}.txt"


def generate_incomplete_cascade(cascade_n: int, removed_inputs: int) -> str:
    """Generate (if missing) the monomer file for cascade-n with first K
    inputs removed. K=1 file already ships; K>=2 builds from the full file."""
    fname = _monomer_filename(cascade_n, removed_inputs)
    path  = os.path.join(EXAMPLE_TBNS_DIR, fname)
    if os.path.exists(path):
        return path

    full_path = os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{cascade_n}.txt")
    with open(full_path) as f:
        lines = f.readlines()

    out_lines = []
    dropped   = 0
    for line in lines:
        if dropped < removed_inputs and line.strip() and not line.startswith("#"):
            dropped += 1
            continue
        out_lines.append(line)
    with open(path, "w") as f:
        f.writelines(out_lines)
    return path


if __name__ == "__main__":
    for removed in [1, 2, 3, 4]:
        p = generate_incomplete_cascade(7, removed)
        eps = build_expected_polymers(7, removed_inputs=removed)
        print(f"\ncascade-7, removed_inputs={removed} -> {os.path.basename(p)}")
        print(f"  |expected| = {len(eps)} polymers")
        for ep in eps[:3]:
            nz = [i for i, c in enumerate(ep.vector) if c]
            print(f"  {ep.label:32s} c={ep.concentration:.1e} M  monomer-idx={nz}")
        print(f"  ... ({len(eps)-4} more)")
        last = eps[-1]
        nz = [i for i, c in enumerate(last.vector) if c]
        print(f"  {last.label:32s} c={last.concentration:.1e} M  monomer-idx={nz}")
