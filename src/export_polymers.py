"""
export_polymers.py  —  Save Pareto-optimal polymer vectors to a flat text file.

After running covering_pipeline.py you have a Python set of tuples in memory.
This helper shows how to write them to the format expected by run_coffee.py,
and can also be imported and called directly from covering_pipeline.py if you
want to save vectors automatically at the end of a run.

Standalone usage (if you serialised the vectors to a .pkl file):
    python export_polymers.py --input vectors.pkl --output pareto_polymers.txt

Import usage (call from covering_pipeline.py):
    from export_polymers import save_polymer_vectors
    save_polymer_vectors(all_hilbert_vectors, "pareto_polymers.txt")

File format produced
--------------------
One line per polymer. Each line is a space-separated list of non-negative
integers of length n_monomers, giving the copy count of each monomer.
Lines where all values are zero are omitted.
Comment lines starting with # are ignored by run_coffee.py.

Example (3 monomers, 4 polymers):
    # Pareto-optimal polymer vectors — n_monomers=3
    1 0 1
    2 1 0
    0 1 2
    1 1 1
"""

import argparse
from curses import raw
import pickle
import sys


def save_polymer_vectors(
    vectors: set[tuple[int, ...]],
    output_path: str,
    n_monomers: int | None = None,
    comment: str = "",
) -> None:
    """
    Parameters
    ----------
    vectors : set of tuples
        The all_hilbert_vectors set from covering_pipeline.py.
    output_path : str
        Destination file path.
    n_monomers : int, optional
        Written into the header comment for documentation. Inferred from
        the first vector if not provided.
    comment : str, optional
        Extra text appended to the header comment line.
    """
    non_zero = [v for v in vectors if any(x != 0 for x in v)]
    if not non_zero:
        print("Warning: no non-zero polymer vectors to write.", file=sys.stderr)
        return

    if n_monomers is None:
        n_monomers = len(non_zero[0])

    bad_lengths = sorted({len(v) for v in non_zero if len(v) != n_monomers})
    if bad_lengths:
        raise ValueError(
            f"save_polymer_vectors length mismatch: "
            f"n_monomers={n_monomers}, bad vector lengths={bad_lengths}, "
            f"output_path={output_path}"
        )

    # Sort for deterministic output (lexicographic by vector)
    non_zero.sort()

    with open(output_path, "w") as f:
        header = f"# Pareto-optimal polymer vectors — n_monomers={n_monomers}"
        if comment:
            header += f" — {comment}"
        f.write(header + "\n")

        for vec in non_zero:
            f.write(" ".join(str(x) for x in vec) + "\n")

    print(f"Saved {len(non_zero)} polymer vectors to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export Pareto-optimal polymer vectors from a pickle file to a flat text file.",
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to a .pkl file containing a set of polymer vector tuples."
    )
    parser.add_argument(
        "--output", default="pareto_polymers.txt",
        help="Output file path. Default: pareto_polymers.txt"
    )
    args = parser.parse_args()

    with open(args.input, "rb") as f:
        raw = pickle.load(f)

    if not isinstance(raw, (set, list)):
        print(f"Error: expected a set or list in {args.input}, got {type(raw)}.", file=sys.stderr)
        sys.exit(1)

    vectors: set[tuple[int, ...]] = set(raw)

    save_polymer_vectors(vectors, args.output)


if __name__ == "__main__":
    main()