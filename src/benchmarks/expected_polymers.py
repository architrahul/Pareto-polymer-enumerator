"""
expected_polymers.py — Build the "expected" polymer set + concentrations
for a linear-cascade system under leakage (one or more inputs removed).

The leakage scenario removes input monomers from the monomer file. The expected
reference used by the leakage experiments is: removed inputs have zero expected
concentration, surviving inputs are free singletons at the initial concentration,
and each module sits in its internal leakage polymers.

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
    removed_input: int | None = None,
) -> list[ExpectedPolymer]:
    """Return expected leakage polymers for cascade-n.

    Backward-compatible mode: ``removed_inputs=K`` means the first K inputs
    are removed and all modules use the unshifted leakage configuration.

    Single-input mode: ``removed_input=j`` means only x_j is removed. The
    expected reference follows the cascade logic:
      * module i is shifted iff both of its inputs are available, which for a
        single removed x_j means i < j-1;
      * shifted module i forms
          (input_a, input_b, m1), (m2, m3), (m4), (m5, m6), output m7;
        if the next module is also shifted, that output m7 is consumed as
        input_a of the next shifted module; otherwise it remains a singleton;
      * unshifted modules form (m1,m2), (m3,m4,m5), (m6,m7).
    """
    total_inputs = cascade_n + 1
    if removed_input is not None:
        if removed_input < 1 or removed_input > total_inputs:
            raise ValueError(f"removed_input must be in [1, {total_inputs}], got {removed_input}")
        removed_set = {removed_input}
        path = generate_single_input_removed_cascade(cascade_n, removed_input)
        shifted_modules = set(range(1, max(1, removed_input - 1)))
    else:
        if removed_inputs < 1 or removed_inputs > total_inputs:
            raise ValueError(f"removed_inputs must be in [1, {total_inputs}], got {removed_inputs}")
        removed_set = set(range(1, removed_inputs + 1))
        path = generate_incomplete_cascade(cascade_n, removed_inputs)
        shifted_modules = set()

    monomers = parse_monomers(path)
    n_total = len(monomers)
    present_inputs = [i for i in range(1, total_inputs + 1) if i not in removed_set]
    input_vec_index = {input_idx: pos for pos, input_idx in enumerate(present_inputs)}
    module_base = len(present_inputs)

    def midx(module_i: int, offset: int) -> int:
        return module_base + (module_i - 1) * MONOMERS_PER_MODULE + offset

    expected_polys: list[ExpectedPolymer] = []
    used = [0] * n_total

    def add_poly(label: str, indices: list[int], conc: float = initial_conc):
        vec = [0] * n_total
        for idx in indices:
            if idx < 0 or idx >= n_total:
                raise RuntimeError(f"expected polymer index {idx} out of range for {label}")
            vec[idx] = 1
            used[idx] += 1
        expected_polys.append(ExpectedPolymer(label=label, vector=tuple(vec), concentration=conc))

    # Surviving external inputs that are not consumed by a shifted module stay
    # as free singleton polymers.
    consumed_external_inputs: set[int] = set()

    for mod_i in range(1, cascade_n + 1):
        if mod_i in shifted_modules:
            # Builder offsets: m2=0, m1=1, m4=2, m5=3, m3=4, m6=5, m7=6.
            if mod_i == 1:
                input_a_idx = input_vec_index[1]
                input_a_label = "x1"
                consumed_external_inputs.add(1)
            else:
                input_a_idx = midx(mod_i - 1, 6)  # previous module output m7
                input_a_label = f"y{mod_i-1}"

            ext_input = mod_i + 1
            input_b_idx = input_vec_index[ext_input]
            consumed_external_inputs.add(ext_input)

            add_poly(f"module{mod_i}_shifted_input_complex",
                     [input_a_idx, input_b_idx, midx(mod_i, 1)])  # i1, i2, m1
            add_poly(f"module{mod_i}_shifted_m2_m3",
                     [midx(mod_i, 0), midx(mod_i, 4)])
            add_poly(f"module{mod_i}_shifted_m4",
                     [midx(mod_i, 2)])
            add_poly(f"module{mod_i}_shifted_m5_m6",
                     [midx(mod_i, 3), midx(mod_i, 5)])

            # Output m7 is consumed by the next shifted module if there is one;
            # otherwise it remains as the free output of this shifted prefix.
            if (mod_i + 1) not in shifted_modules:
                add_poly(f"module{mod_i}_shifted_output", [midx(mod_i, 6)])
        else:
            for poly_idx, offsets in enumerate(MODULE_POLYMER_OFFSETS):
                kind = ["central_complement", "intermediate", "output_pair"][poly_idx]
                add_poly(f"module{mod_i}_{kind}", [midx(mod_i, off) for off in offsets])

    for input_idx in present_inputs:
        if input_idx in consumed_external_inputs:
            continue
        add_poly(f"input_x{input_idx}", [input_vec_index[input_idx]])

    missing = [i for i, c in enumerate(used) if c == 0]
    duplicated = [i for i, c in enumerate(used) if c > 1]
    if missing or duplicated:
        raise RuntimeError(
            f"expected polymer construction invalid: missing={missing[:8]}, "
            f"duplicated={duplicated[:8]} (removed_input={removed_input}, removed_inputs={removed_inputs})"
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


def _single_removed_monomer_filename(cascade_n: int, removed_input: int) -> str:
    if removed_input == 1:
        return f"monomers_cascade_n{cascade_n}_incomplete.txt"
    return f"monomers_cascade_n{cascade_n}_missing_input{removed_input}.txt"


def generate_single_input_removed_cascade(cascade_n: int, removed_input: int) -> str:
    """Generate/read cascade-n with exactly one input x_removed_input removed."""
    total_inputs = cascade_n + 1
    if removed_input < 1 or removed_input > total_inputs:
        raise ValueError(f"removed_input must be in [1, {total_inputs}], got {removed_input}")

    fname = _single_removed_monomer_filename(cascade_n, removed_input)
    path = os.path.join(EXAMPLE_TBNS_DIR, fname)
    if os.path.exists(path):
        return path

    full_path = os.path.join(EXAMPLE_TBNS_DIR, f"monomers_cascade_n{cascade_n}.txt")
    with open(full_path) as f:
        lines = f.readlines()

    out_lines = []
    input_seen = 0
    for line in lines:
        stripped = line.strip()
        if stripped and not line.startswith("#") and input_seen < total_inputs:
            input_seen += 1
            if input_seen == removed_input:
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
