import subprocess
import time
from datetime import datetime
import os
import sys
import argparse
import itertools
import math
import random
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import threading
import atexit
import glob
import json
import sqlite3
from functools import lru_cache
from export_polymers import save_polymer_vectors
from paths import DEFAULT_MONOMER_FILE, NORMALIZ_EXE, LOGS_DIR, RESULTS_DIR, COVERING_DP_DB

# -------------------------
# Config
# -------------------------

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PYTHON_SCRIPT = os.path.join(SRC_DIR, "monomers_to_normaliz.py")

SCRATCH_DIR  = os.path.join(SRC_DIR, ".normaliz_tmp")
TMP_MONOMERS = os.path.join(SCRATCH_DIR, "tmp_monomers.txt")
TMP_VECTORS  = os.path.join(SCRATCH_DIR, "vectors.txt")
TMP_EQS      = os.path.join(SCRATCH_DIR, "eqs")
PROBE_LIMIT  = 100
K_MAX        = 25

# Any single block exceeding this is killed and treated as contributing zero vectors to the union 
NORMALIZ_TIMEOUT_SECONDS = 3 * 3600  # 3 hours

# -------------------------
# Utility functions
# -------------------------

def ensure_scratch_dir():
    """Create the shared Normaliz scratch directory under src/."""
    os.makedirs(SCRATCH_DIR, exist_ok=True)


def cleanup_normaliz_files(remove_dir=False):
    """Remove generated Normaliz/converter scratch files.

    This repo assumes one Normaliz run at a time. Scratch files live directly in
    src/.normaliz_tmp so even interrupted runs leave only a small bounded set of
    files, which the next block removes before starting.
    """
    if not os.path.isdir(SCRATCH_DIR):
        return
    for path in glob.glob(os.path.join(SCRATCH_DIR, "eqs.*")):
        if os.path.isfile(path):
            os.remove(path)
    for name in ["tmp_monomers.txt", "vectors.txt"]:
        path = os.path.join(SCRATCH_DIR, name)
        if os.path.exists(path):
            os.remove(path)


atexit.register(cleanup_normaliz_files)


def read_hilbert_basis(base_filename=TMP_EQS):
    output_file = f"{base_filename}.out"
    if not os.path.exists(output_file):
        return []

    vectors = []
    with open(output_file, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "Hilbert basis elements:" in line:
            start = i + 1
            break
    if start is None:
        return []

    for line in lines[start:]:
        line = line.strip()
        if not line or line.startswith("***") or "extreme rays" in line.lower():
            break
        if line.startswith("#"):
            continue
        try:
            vectors.append([int(x) for x in line.split()])
        except ValueError:
            pass

    return vectors


def expand_vector_to_full_space(reduced_vector, selected_indices, n):
    """Monomer mode: maps reduced monomer-subset vector back to full monomer space."""
    full_vector = [0] * n
    for i, idx in enumerate(selected_indices):
        if i < len(reduced_vector):
            full_vector[idx] = reduced_vector[i]
    return tuple(full_vector)


def expand_vector_to_full_monomer_space(reduced_vector, filtered_monomers_indices, n_monomers):
    """Domain mode: maps filtered-monomer vector back to full monomer space."""
    monomer_part = reduced_vector[:len(filtered_monomers_indices)]
    full_vector  = [0] * n_monomers
    for i, idx in enumerate(filtered_monomers_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    return tuple(full_vector)


_skip_current = False

def _listen_for_skip():
    global _skip_current
    while True:
        try:
            if input().strip().lower() == 's':
                _skip_current = True
        except EOFError:
            # No interactive stdin (running headless under nohup/caffeinate).
            # No skip key available; exit the listener silently.
            return

def start_input_listener():
    threading.Thread(target=_listen_for_skip, daemon=True).start()

def check_and_clear_skip() -> bool:
    global _skip_current
    if _skip_current:
        _skip_current = False
        return True
    return False


# -------------------------
# Domain helpers
# -------------------------

def get_domains_from_monomer(monomer: str) -> set:
    return {domain.rstrip("*") for domain in monomer.split()}


def get_all_unique_domains(monomers: list) -> list:
    seen = OrderedDict()
    for monomer in monomers:
        for domain in get_domains_from_monomer(monomer):
            seen[domain] = None
    return list(seen.keys())


def filter_monomers_by_domains(monomers: list, selected_domains: list) -> list:
    """Return only monomers whose domains are all within the selected set."""
    selected_set = set(selected_domains)
    return [m for m in monomers if get_domains_from_monomer(m).issubset(selected_set)]


# -------------------------
# Covering design fetch / compute
# -------------------------

_LJCR_SESSION: "requests.Session | None" = None
_LJCR_SESSION_LOCK = threading.Lock()


def _ljcr_session():
    """Lazily build a requests.Session so multiple fetches share one keep-alive
    TCP connection — turns ~300ms/fetch into ~50ms/fetch after the first."""
    global _LJCR_SESSION
    if _LJCR_SESSION is None:
        with _LJCR_SESSION_LOCK:
            if _LJCR_SESSION is None:
                s = requests.Session()
                s.headers.update({"User-Agent": "Mozilla/5.0"})
                _LJCR_SESSION = s
    return _LJCR_SESSION


def fetch_covering_online(v: int, k: int, t: int) -> list:
    url = f"https://ljcr.dmgordon.org/show_cover.php?v={v}&k={k}&t={t}"
    r   = _ljcr_session().get(url, timeout=30)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "lxml")
    pre  = soup.find("pre")
    if pre is None:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")

    blocks = []
    for line in pre.text.strip().splitlines():
        row = [int(x) for x in line.split()]
        if row:
            blocks.append(row)

    if not blocks:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")
    return blocks


# SQLite recipe support for --fallback-dp.
#
# The bundled DB stores compact GPK Section 5 construction choices. It does not
# store all covering blocks; blocks are materialized on demand, fetching seed
# designs from LJCR when needed.
#   COVERING_DP_DB          optional override for the recipe DB
#   COVERING_REPO_DB        optional local seed-block cache
#   COVERING_DP_RUN_NAME    optional run_name
_OFFLINE_COVER_CACHE: dict = {}


def _normalize_blocks(blocks):
    seen = set()
    out = []
    for block in blocks:
        b = tuple(sorted(int(x) for x in block))
        if b not in seen:
            seen.add(b)
            out.append(list(b))
    return out


def _covering_recipe_paths():
    # Default to the bundled compact DP recipe DB. Users can override with
    # COVERING_DP_DB after rebuilding the table. COVERING_REPO_DB is optional;
    # without it, seed designs are fetched from LJCR online as needed.
    dp_db = os.environ.get("COVERING_DP_DB") or COVERING_DP_DB
    if not dp_db or not os.path.exists(dp_db):
        return None
    repo_db = os.environ.get("COVERING_REPO_DB")
    if repo_db and not os.path.exists(repo_db):
        repo_db = None
    run_name = os.environ.get("COVERING_DP_RUN_NAME") or "gpk_V150_K80_T8"
    return repo_db, dp_db, run_name


def _trivial_cover(key):
    v, k, t = key
    if t == 0:
        return [list(range(1, k + 1))]
    if k == v:
        return [list(range(1, v + 1))]
    if t == 1:
        blocks = []
        for start in range(1, v + 1, k):
            core = list(range(start, min(start + k, v + 1)))
            filler = [x for x in range(1, v + 1) if x not in core]
            blocks.append(sorted(core + filler[: k - len(core)]))
        return _normalize_blocks(blocks)
    if t == k:
        return [list(c) for c in itertools.combinations(range(1, v + 1), k)]
    raise RuntimeError(f"No trivial covering construction for C{key}")


def _fetch_seed_blocks(repo_conn, key):
    """Fetch seed blocks from an optional local repo DB, else LJCR online."""
    if repo_conn is not None:
        rows = repo_conn.execute(
            "SELECT block FROM seed_design_blocks WHERE v=? AND k=? AND t=? ORDER BY block_index",
            key,
        ).fetchall()
        if rows:
            return [[int(x) for x in row[0].split()] for row in rows]
    return fetch_covering_online(*key)


def _recipe_run_names(dp_conn, run_name, key):
    if run_name:
        return [run_name]
    rows = dp_conn.execute(
        """
        SELECT b.run_name
        FROM dp_bounds b
        LEFT JOIN dp_runs r ON r.run_name = b.run_name
        WHERE b.v=? AND b.k=? AND b.t=?
        ORDER BY COALESCE(r.created_at, '') DESC, b.run_name DESC
        """,
        key,
    ).fetchall()
    return [row[0] for row in rows]


def _try_covering_dp_recipe(v: int, k: int, t: int):
    """Materialize blocks from the bundled GPK DP recipe DB, or None."""
    key = (v, k, t)
    if key in _OFFLINE_COVER_CACHE:
        return _OFFLINE_COVER_CACHE[key]

    paths = _covering_recipe_paths()
    if paths is None:
        return None
    repo_db, dp_db, requested_run = paths

    try:
        repo_conn = sqlite3.connect(repo_db) if repo_db else None
        dp_conn = sqlite3.connect(dp_db)
        if repo_conn is not None:
            repo_conn.row_factory = sqlite3.Row
        dp_conn.row_factory = sqlite3.Row
    except sqlite3.Error:
        return None

    try:
        run_names = _recipe_run_names(dp_conn, requested_run, key)
        if not run_names:
            return None

        def product(left_blocks, right_blocks, right_offset, local_v, local_k):
            out = []
            for left in left_blocks:
                for right in right_blocks:
                    merged = sorted(set(left) | {x + right_offset for x in right})
                    if len(merged) < local_k:
                        used = set(merged)
                        filler = [x for x in range(1, local_v + 1) if x not in used]
                        merged = sorted(merged + filler[: local_k - len(merged)])
                    out.append(merged)
            return _normalize_blocks(out)

        for run_name in run_names:
            @lru_cache(maxsize=None)
            def build_cover(cur_key):
                row = dp_conn.execute(
                    "SELECT * FROM dp_bounds WHERE run_name=? AND v=? AND k=? AND t=?",
                    (run_name, cur_key[0], cur_key[1], cur_key[2]),
                ).fetchone()
                if row is None:
                    # A direct seed may exist even if no reason row was stored.
                    return tuple(tuple(b) for b in _fetch_seed_blocks(repo_conn, cur_key))

                rtype = row["reason_kind"]
                if rtype == "seed":
                    if row["source"] == "trivial":
                        return tuple(tuple(b) for b in _trivial_cover(cur_key))
                    return tuple(tuple(b) for b in _fetch_seed_blocks(repo_conn, cur_key))
                if rtype == "lift":
                    src = (row["src_v"], row["src_k"], row["src_t"])
                    lifted = []
                    for block in build_cover(src):
                        b = sorted(set(block) | {cur_key[0]})
                        if len(b) < cur_key[1]:
                            used = set(b)
                            filler = [x for x in range(1, cur_key[0] + 1) if x not in used]
                            b = sorted(b + filler[: cur_key[1] - len(b)])
                        lifted.append(b)
                    return tuple(tuple(b) for b in lifted)
                if rtype == "weaken":
                    src = (row["src_v"], row["src_k"], row["src_t"])
                    return build_cover(src)
                if rtype == "section5":
                    tree = json.loads(row["tree_json"])

                    def build_interval(node):
                        if node["type"] == "direct":
                            left_key = tuple(node["left"])
                            right_key = tuple(node["right"])
                            return product(
                                [list(b) for b in build_cover(left_key)],
                                [list(b) for b in build_cover(right_key)],
                                left_key[0],
                                cur_key[0],
                                cur_key[1],
                            )
                        return _normalize_blocks(build_interval(node["left"]) + build_interval(node["right"]))

                    return tuple(tuple(b) for b in build_interval(tree))
                raise RuntimeError(f"Unknown fallback-DP recipe reason kind {rtype!r}")

            try:
                blocks = [list(b) for b in build_cover(key)]
                if blocks:
                    _OFFLINE_COVER_CACHE[key] = blocks
                    print(f"  Built fallback-DP covering C({v},{k},{t}) with {len(blocks)} blocks.")
                    return blocks
            except Exception as e:
                raise RuntimeError(
                    f"Failed to materialize C({v},{k},{t}) from covering DP database. "
                    "The fallback construction requires internet access to fetch seed designs from LJCR."
                ) from e
    finally:
        try:
            if repo_conn is not None:
                repo_conn.close()
            dp_conn.close()
        except Exception:
            pass
    return None


# In-memory caches (per process; never persisted):
#   _DP_SIZE  : (v, k, t) -> integer size of the design _DP would produce
#   _DP_PLAN  : (v, k, t) -> 'plan' explaining how to build the design (see below)
#   _LJCR_CACHE : (v, k, t) -> actual block list fetched from LJCR
#   _COVERING_DP_CACHE : (v, k, t) -> materialised block list (final designs)
#
# Two-pass strategy (mirrors Minki Hhan's vibecoded approach + GPK §5):
#   pass 1: _dp_size() recursively computes upper-bound sizes via the
#           recurrence, picking the best split parameters. No blocks are
#           constructed yet; LJCR is hit at most once per (v,k,t) leaf.
#   pass 2: compute_covering_dp() walks the saved plans to materialise
#           actual blocks on the chosen path only.
_DP_SIZE: dict = {}
_DP_PLAN: dict = {}
_LJCR_CACHE: dict = {}
_COVERING_DP_CACHE: dict = {}

# (v // 2, ±1, ±2, plus v/3 and 2v/3) — a small balanced set of splits.
def _split_candidates(v: int, t: int):
    raw = {v // 2 - 2, v // 2 - 1, v // 2, v // 2 + 1, v // 2 + 2,
           v // 3, 2 * v // 3, t, v - t}
    return sorted(x for x in raw if 1 <= x < v)


def _is_ljcr_leaf(v: int, k: int, t: int) -> bool:
    return k <= 25 and t <= 8 and 1 <= v < 100


def _prefetch_ljcr_concurrent(triples, max_workers: int = 32):
    """Concurrently fetch LJCR designs for any (v,k,t) we haven't yet seen.
    Successes go into _LJCR_CACHE; failures get a None sentinel so we don't
    retry. Populates _DP_SIZE / _DP_PLAN for any that succeed.
    """
    todo = [tr for tr in set(triples)
            if _is_ljcr_leaf(*tr) and tr not in _LJCR_CACHE]
    if not todo:
        return

    def _one(tr):
        try:
            return tr, fetch_covering_online(*tr)
        except (RuntimeError, requests.RequestException):
            return tr, None

    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        for tr, blocks in pool.map(_one, todo):
            _LJCR_CACHE[tr] = blocks
            if blocks is not None:
                _DP_SIZE[tr] = len(blocks)
                _DP_PLAN[tr] = ("ljcr",)


def _enumerate_subqueries(v: int, k: int, t: int):
    """List of (v', k', t') children of _dp_size(v, k, t)'s refined-DP
    iteration — used to know what to prefetch before recursing.

    The refined construction asks for C(v1, ℓ, j) · C(v2, k-ℓ, t-i) for every
    interval [i, j] ⊆ [0, t]. The j-values run over {0..t}; the t-i values
    also run over {0..t} (since 0 ≤ i ≤ j ≤ t).
    """
    out = []
    for v1 in _split_candidates(v, t):
        v2 = v - v1
        for j in range(t + 1):
            for i in range(j + 1):
                ell_lo = max(j, k - v2)
                ell_hi = min(v1, k - (t - i))
                for ell in range(ell_lo, ell_hi + 1):
                    out.append((v1, ell, j))
                    out.append((v2, k - ell, t - i))
    return out


def _dp_size(v: int, k: int, t: int) -> int:
    """Return upper bound on C(v, k, t) and memoise the best construction plan
    in _DP_PLAN. Returns math.inf when no construction is possible.

    Plans are tagged tuples consumed by `_materialise`:
      ("trivial", "all")     -> single full v-element block
      ("trivial", "comb")    -> every t-subset of [1..v]
      ("trivial", "pack")    -> ⌈v/k⌉ disjoint k-windows (t = 1)
      ("trivial", "empty")   -> one empty block (t = 0, k = 0)
      ("trivial", "one")     -> one [1..k] block (t = 0, k > 0)
      ("ljcr",)              -> blocks cached in _LJCR_CACHE
      ("weaken", src)        -> reuse a C(v,k,t+1) cover as a C(v,k,t) cover
      ("lift", src)          -> append point v to every C(v-1,k-1,t) block
      ("split", v1, root)    -> GPK §5 interval-DP split
    """
    if (v, k, t) in _DP_SIZE:
        return _DP_SIZE[(v, k, t)]

    # --- Trivial base cases ---------------------------------------------
    if t == 0:
        plan = ("trivial", "empty" if k <= 0 else "one")
        _DP_SIZE[(v, k, t)] = 1
        _DP_PLAN[(v, k, t)] = plan
        return 1

    if k <= 0 or k < t:
        _DP_SIZE[(v, k, t)] = math.inf
        return math.inf

    if k >= v:
        _DP_SIZE[(v, k, t)] = 1
        _DP_PLAN[(v, k, t)] = ("trivial", "all")
        return 1

    if t == 1:
        size = (v + k - 1) // k
        _DP_SIZE[(v, k, t)] = size
        _DP_PLAN[(v, k, t)] = ("trivial", "pack")
        return size

    if k == t:
        size = math.comb(v, t)
        _DP_SIZE[(v, k, t)] = size
        _DP_PLAN[(v, k, t)] = ("trivial", "comb")
        return size

    # --- LJCR leaf ------------------------------------------------------
    if _is_ljcr_leaf(v, k, t):
        if (v, k, t) in _LJCR_CACHE:
            cached = _LJCR_CACHE[(v, k, t)]
            if cached is not None:
                _DP_SIZE[(v, k, t)] = len(cached)
                _DP_PLAN[(v, k, t)] = ("ljcr",)
                return len(cached)
        else:
            try:
                blocks = fetch_covering_online(v, k, t)
                _LJCR_CACHE[(v, k, t)] = blocks
                _DP_SIZE[(v, k, t)] = len(blocks)
                _DP_PLAN[(v, k, t)] = ("ljcr",)
                return len(blocks)
            except (RuntimeError, requests.RequestException):
                _LJCR_CACHE[(v, k, t)] = None   # cache the "miss"

    # --- Recursive improvements used by the fallback-DP recipe code -----
    # Weaken: a (v,k,t+1)-cover also covers every t-subset.
    # Lift: from C(v-1,k-1,t), append the new point v to every block.
    # These two rules are especially important for k > 25 because they can
    # reduce quickly to an LJCR seed instead of invoking an expensive Section 5
    # construction from scratch.
    best = math.inf
    best_plan = None

    if t + 1 <= k:
        src = (v, k, t + 1)
        cand = _dp_size(*src)
        if cand < best:
            best = cand
            best_plan = ("weaken", src)

    if v >= 1 and k >= 1:
        src = (v - 1, k - 1, t)
        cand = _dp_size(*src)
        if cand < best:
            best = cand
            best_plan = ("lift", src)

    # --- GPK §5 refined c_{i,j} interval DP -----------------------------
    # Prefetch all leaf sub-queries concurrently before iterating.
    subs = _enumerate_subqueries(v, k, t)
    _prefetch_ljcr_concurrent(subs)
    for v1 in _split_candidates(v, t):
        v2 = v - v1

        # c[(i,j)] = (size, node) for 0 ≤ i ≤ j ≤ t.
        # node = ("prod", j, t-i, ℓ)                — one Cartesian product
        #      | ("split_int", r, left_node, right_node) — split at r ∈ [i, j-1]
        c: dict = {}
        feasible = True
        for j in range(t + 1):
            if not feasible:
                break
            for i in range(j, -1, -1):
                # Option A: one product covering the whole interval [i, j].
                node = None
                size_here = math.inf
                ell_lo = max(j, k - v2)
                ell_hi = min(v1, k - (t - i))
                for ell in range(ell_lo, ell_hi + 1):
                    a = _dp_size(v1, ell, j)
                    b = _dp_size(v2, k - ell, t - i)
                    if a == math.inf or b == math.inf:
                        continue
                    term = a * b
                    if term < size_here:
                        size_here = term
                        node = ("prod", j, t - i, ell)

                # Option B: split [i, j] at r, sum the two sub-intervals.
                for r in range(i, j):
                    left  = c.get((i,     r))
                    right = c.get((r + 1, j))
                    if left is None or right is None:
                        continue
                    if left[0] == math.inf or right[0] == math.inf:
                        continue
                    cand = left[0] + right[0]
                    if cand < size_here:
                        size_here = cand
                        node = ("split_int", r, left[1], right[1])

                if node is None:
                    feasible = False
                    break
                c[(i, j)] = (size_here, node)

        if not feasible or (0, t) not in c:
            continue
        total_size, root_node = c[(0, t)]
        if total_size < best:
            best = total_size
            best_plan = ("split", v1, root_node)

    _DP_SIZE[(v, k, t)] = best
    if best_plan is not None:
        _DP_PLAN[(v, k, t)] = best_plan
    return best


def _materialise(v: int, k: int, t: int) -> list:
    """Build the actual block list from the saved plan."""
    if (v, k, t) in _COVERING_DP_CACHE:
        return _COVERING_DP_CACHE[(v, k, t)]

    plan = _DP_PLAN.get((v, k, t))
    if plan is None:
        raise RuntimeError(f"no plan recorded for C({v}, {k}, {t}); did _dp_size succeed?")

    kind = plan[0]
    if kind == "trivial":
        subkind = plan[1]
        if subkind == "empty":
            result = [[]]
        elif subkind == "one":
            result = [list(range(1, k + 1))]
        elif subkind == "all":
            result = [list(range(1, v + 1))]
        elif subkind == "pack":
            n_blocks = (v + k - 1) // k
            result = []
            for i in range(n_blocks):
                lo = min(i * k, v - k) + 1
                result.append(list(range(lo, lo + k)))
        elif subkind == "comb":
            result = [list(c) for c in itertools.combinations(range(1, v + 1), t)]
        else:
            raise RuntimeError(f"unknown trivial subkind {subkind}")
    elif kind == "ljcr":
        result = _LJCR_CACHE[(v, k, t)]
    elif kind == "weaken":
        _, src = plan
        result = [list(block) for block in _materialise(*src)]
    elif kind == "lift":
        _, src = plan
        result = []
        for block in _materialise(*src):
            lifted = sorted(set(block) | {v})
            if len(lifted) < k:
                used = set(lifted)
                filler = [x for x in range(1, v + 1) if x not in used]
                lifted = sorted(lifted + filler[: k - len(lifted)])
            result.append(lifted)
    elif kind == "split":
        _, v1, root = plan
        v2 = v - v1

        def _walk(node):
            """Materialise blocks for one interval-DP plan node."""
            if node[0] == "prod":
                _, j_val, t_minus_i, ell = node
                D1 = _materialise(v1, ell, j_val)
                D2 = _materialise(v2, k - ell, t_minus_i)
                blocks = []
                for B1 in D1:
                    for B2 in D2:
                        merged = sorted(set(B1) | {x + v1 for x in B2})
                        if len(merged) < k:
                            used   = set(merged)
                            unused = [x for x in range(1, v + 1) if x not in used]
                            merged = sorted(merged + unused[: k - len(merged)])
                        blocks.append(merged)
                return blocks
            if node[0] == "split_int":
                _, _r, left, right = node
                return _walk(left) + _walk(right)
            raise RuntimeError(f"unknown plan node {node[0]}")

        result = _walk(root)
    else:
        raise RuntimeError(f"unknown plan kind {kind}")

    _COVERING_DP_CACHE[(v, k, t)] = result
    return result


def compute_covering_dp(v: int, k: int, t: int) -> list:
    """Build C(v,k,t) from the bundled GPK recipe database.

    This is intentionally strict: if the recipe DB is missing, the requested
    state is outside its range, or an online LJCR seed fetch fails, raise a
    clear error. Do not silently fall back to the older runtime Section-5 search,
    because that can generate different covering sizes and change experiments.
    """
    if (v, k, t) in _COVERING_DP_CACHE:
        return _COVERING_DP_CACHE[(v, k, t)]

    blocks = _try_covering_dp_recipe(v, k, t)
    if blocks is None:
        raise RuntimeError(
            f"Could not build C({v},{k},{t}) with --fallback-dp. "
            "Check that data/covering_design/gpk_dp.sqlite exists, the requested "
            "parameters are covered by it, and internet access to LJCR is available."
        )

    _COVERING_DP_CACHE[(v, k, t)] = blocks
    _DP_SIZE[(v, k, t)] = len(blocks)
    _DP_PLAN[(v, k, t)] = ("fallback_dp_recipe",)
    return blocks


def load_covering_blocks(v: int, k: int, t: int, fallback_dp: bool = False) -> list:
    """
    Returns list of blocks for C(v, k, t).
    Raises RuntimeError if unavailable and fallback_dp is False.
    When k == v the single full block is returned directly.

    If fallback_dp is True and LJCR has no covering for (v, k, t), we
    build one locally via the GPK dynamic-programming construction
    (`compute_covering_dp`).
    """
    if k == v:
        return [list(range(1, v + 1))]

    try:
        return fetch_covering_online(v, k, t)
    except (RuntimeError, requests.RequestException) as e:
        if not fallback_dp:
            raise RuntimeError(
                f"Online fetch failed for C({v},{k},{t}): {e}. "
                f"Use --fallback-dp to build the covering locally."
            ) from e
        print(f"  Online fetch failed ({e}). Trying --fallback-dp construction.")
        return compute_covering_dp(v, k, t)


# -------------------------
# Single-block Normaliz runner
# -------------------------

def run_normaliz_on_subset(subset_monomers: list) -> tuple:
    """
    Write monomers, invoke monomers_to_normaliz.py, run Normaliz.
    Returns (elapsed_seconds, hilbert_vectors_raw).

    If Normaliz exceeds NORMALIZ_TIMEOUT_SECONDS, the child is killed and
    this returns (elapsed, []) so the caller treats the block as
    contributing no Hilbert basis vectors.
    """
    cleanup_normaliz_files()
    ensure_scratch_dir()

    with open(TMP_MONOMERS, "w") as f:
        for m in subset_monomers:
            f.write(m + "\n")

    subprocess.run(
        [
            sys.executable, PYTHON_SCRIPT,
            "--input", TMP_MONOMERS,
            "--vectors", TMP_VECTORS,
            "--eqs-in", f"{TMP_EQS}.in",
        ],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
        cwd=SRC_DIR,
    )

    t0 = time.time()
    try:
        subprocess.run(
            [NORMALIZ_EXE, "-d", "-N", TMP_EQS],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=NORMALIZ_TIMEOUT_SECONDS,
            cwd=SRC_DIR,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        print(
            f"  [WARN] Normaliz exceeded {NORMALIZ_TIMEOUT_SECONDS}s "
            f"({elapsed:.1f}s elapsed); block killed, treated as 0 vectors."
        )
        return elapsed, []
    return time.time() - t0, read_hilbert_basis()


# -------------------------
# Probe phase for one k value
# -------------------------

def probe_k(k, t, blocks, n, all_monomers, mode, all_domains, n_monomers,
            best_projected, log):
    """
    Run up to PROBE_LIMIT randomly sampled blocks for the given k.

    Early-abort rule: if cumulative probe wall time exceeds best_projected at
    any point during the probe, this k is pruned immediately.

    Returns:
        projected_total  – avg_block_time * num_blocks, or None if pruned/empty
        probe_times      – list of individual block times recorded
        num_blocks       – total number of blocks in the covering design
    """
    num_blocks  = len(blocks)
    probe_size  = min(PROBE_LIMIT, num_blocks)
    sample_idxs = random.sample(range(num_blocks), probe_size)

    probe_times = []
    log.write(f"    Probing k={k}: {num_blocks} blocks total, sampling {probe_size}\n")
    log.flush()

    for i, idx in enumerate(sample_idxs):
        if check_and_clear_skip():
            log.write(f"    k={k}: SKIPPED by user during probe\n")
            log.flush()
            return None, probe_times, num_blocks

        block = blocks[idx]
        cleanup_normaliz_files()

        if mode == "monomer":
            selected_indices = [x - 1 for x in block]
            subset           = [all_monomers[j] for j in selected_indices]
        else:
            selected_domains = [all_domains[x - 1] for x in block]
            subset           = filter_monomers_by_domains(all_monomers, selected_domains)
            if not subset:
                continue

        elapsed, _ = run_normaliz_on_subset(subset)
        probe_times.append(elapsed)
        cumulative = sum(probe_times)

        # Early abort: cumulative probe time already exceeds best projected total
        if best_projected is not None and cumulative > best_projected:
            log.write(
                f"    k={k}: PROBE PRUNED at iteration {i+1}/{probe_size} "
                f"(cumulative={cumulative:.3f}s > best_projected={best_projected:.3f}s)\n"
            )
            log.flush()
            print(
                f"  [t={t}] k={k}: pruned at block {i+1}/{probe_size} "
                f"(cumulative={cumulative:.3f}s > best={best_projected:.3f}s)"
            )
            return None, probe_times, num_blocks

    if not probe_times:
        return None, probe_times, num_blocks

    avg_time        = sum(probe_times) / len(probe_times)
    projected_total = avg_time * num_blocks

    log.write(
        f"    k={k}: probe complete — "
        f"probe_iters={len(probe_times)}, avg={avg_time:.4f}s, "
        f"num_blocks={num_blocks}, projected_total={projected_total:.3f}s\n"
    )
    log.flush()
    return projected_total, probe_times, num_blocks


# -------------------------
# Full enumeration for one k value
# -------------------------

def full_run_k(k, blocks, n, all_monomers, mode, all_domains, n_monomers, log):
    """
    Run Normaliz on every block for the given k, collecting all Hilbert vectors.
    Returns (result_dict, total_normaliz_time).
    """
    num_blocks          = len(blocks)
    all_hilbert_vectors = set()
    times               = []
    wall_start          = time.time()

    log.write(f"    Full run k={k}: {num_blocks} blocks\n")
    log.flush()

    for idx, block in enumerate(blocks):
        if check_and_clear_skip():
            log.write(f"    k={k}: ABORTED by user during full run\n")
            log.flush()
            return None, float("inf")

        cleanup_normaliz_files()

        if mode == "monomer":
            selected_indices = [x - 1 for x in block]
            subset           = [all_monomers[j] for j in selected_indices]
        else:
            selected_domains = [all_domains[x - 1] for x in block]
            subset           = filter_monomers_by_domains(all_monomers, selected_domains)
            if not subset:
                continue
            selected_indices = [
                i for i, m in enumerate(all_monomers) if m in subset
            ]

        elapsed, raw_vectors = run_normaliz_on_subset(subset)
        times.append(elapsed)

        for rv in raw_vectors:
            if mode == "monomer":
                all_hilbert_vectors.add(
                    expand_vector_to_full_space(rv, selected_indices, n)
                )
            else:
                all_hilbert_vectors.add(
                    expand_vector_to_full_monomer_space(rv, selected_indices, n_monomers)
                )

        if (idx + 1) % 50 == 0 or (idx + 1) == num_blocks:
            print(f"  k={k}: {idx+1}/{num_blocks} blocks done "
                  f"({sum(times):.1f}s, {len(all_hilbert_vectors)} vectors)")

    wall_time     = time.time() - wall_start
    normaliz_time = sum(times)
    overhead      = wall_time - normaliz_time
    avg_normaliz  = normaliz_time / num_blocks if num_blocks else 0.0

    log.write(
        f"    k={k}: FULL RUN COMPLETE — "
        f"wall={wall_time:.3f}s  normaliz={normaliz_time:.3f}s  "
        f"overhead={overhead:.3f}s  avg/block={avg_normaliz:.4f}s  "
        f"unique_vectors={len(all_hilbert_vectors)}\n"
    )
    log.flush()

    return {
        "k":                           k,
        "num_blocks":                  num_blocks,
        "total_wall_time":             wall_time,
        "total_normaliz_time":         normaliz_time,
        "overhead_time":               overhead,
        "avg_normaliz_time_per_block": avg_normaliz,
        "unique_vectors":              len(all_hilbert_vectors),
        "vectors":                     all_hilbert_vectors,
    }, normaliz_time


## -------------------------
# Main covering sweep for one t value
# -------------------------

def run_covering_sweep(
    t,
    all_monomers,
    mode,
    log,
    fallback_dp = False,
    probe_only      = False,
    include_base    = False,
    save            = False,
    save_dir        = None,
):
    """
    Sweep k from t+1 up to min(K_MAX, n), probing each k with probe-and-prune.
    Skips any k where no covering design is available (no crash).

    probe_only=True  → returns (best_k, best_projected, None)
    probe_only=False → runs full enumeration on best k;
                       returns (best_k, best_projected, full_result_dict)

    include_base=True additionally probes/runs k=n.
    If include_base=True, k=n is always run FIRST to establish the initial upper bound.
    """
    n_monomers  = len(all_monomers)
    all_domains = get_all_unique_domains(all_monomers)
    n_domains   = len(all_domains)
    n           = n_monomers if mode == "monomer" else n_domains

    # Candidate k values excluding base case
    k_values = list(range(t + 1, min(K_MAX, n) + 1))
    if n in k_values:
        k_values.remove(n)

    log.write(
        f"\n{'='*60}\n"
        f"  Covering sweep: t={t}, mode={mode}, n={n}\n"
        f"  include_base={include_base}\n"
        f"  k_values={([n] if include_base else []) + k_values}\n"
        f"{'='*60}\n"
    )
    log.flush()

    best_projected = None
    best_k         = None
    best_blocks    = None
    probe_summary  = {}  # k -> projected_total or None
    consecutive_increases = 0
    last_projected = None

    # ------------------------------------------------------------
    # STEP 1: Run base case FIRST if requested
    # ------------------------------------------------------------
    if include_base:
        print(f"\n  [t={t}] Probing base case k={n} first ...")

        try:
            base_blocks = load_covering_blocks(n, n, t, fallback_dp=fallback_dp)
        except RuntimeError as e:
            print(f"  Base case k={n} unavailable: {e}")
            log.write(f"    k={n}: SKIPPED (base case unavailable: {e})\n")
            log.flush()
            probe_summary[n] = None
        else:
            projected, probe_times, num_blocks = probe_k(
                n, t, base_blocks, n, all_monomers, mode, all_domains, n_monomers,
                None, log
            )
            probe_summary[n] = projected

            if projected is not None:
                best_projected = projected
                best_k         = n
                best_blocks    = base_blocks
                print(f"  [t={t}] Base case k={n}: projected={projected:.3f}s")

                # --- NEW: full run + save for base case ---
                if not probe_only:
                    print(f"\n  [t={t}] Running full enumeration on base case k={n} ...")
                    base_result, _ = full_run_k(
                        n, base_blocks, n, all_monomers, mode, all_domains, n_monomers, log
                    )
                    if base_result is not None and save and save_dir:
                        os.makedirs(save_dir, exist_ok=True)
                        output_path = os.path.join(
                            save_dir, f"hilbert_k{n}_t{t}_{mode}.txt"
                        )
                        save_polymer_vectors(
                            base_result["vectors"],
                            output_path,
                            n_monomers=n_monomers,
                            comment=f"covering mode, k={n}, t={t}, mode={mode} (base case)"
                        )
                        print(f"  Saved base case vectors → {output_path}")

    # ------------------------------------------------------------
    # STEP 2: Probe remaining k values using current best for pruning
    # ------------------------------------------------------------
    for k in k_values:
        print(f"\n  [t={t}] Probing k={k} ...")

        try:
            blocks = load_covering_blocks(n, k, t, fallback_dp=fallback_dp)
        except RuntimeError as e:
            print(f"  Skipping k={k}: {e}")
            log.write(f"    k={k}: SKIPPED (covering unavailable: {e})\n")
            log.flush()
            probe_summary[k] = None
            continue

        projected, probe_times, num_blocks = probe_k(
            k, t, blocks, n, all_monomers, mode, all_domains, n_monomers,
            best_projected, log
        )

        probe_summary[k] = projected

        if projected is None:
            continue

        if best_projected is None:
            print(f"  [t={t}] k={k}: projected={projected:.3f}s  (first estimate)")
        else:
            print(f"  [t={t}] k={k}: projected={projected:.3f}s  "
                  f"(best so far: {best_projected:.3f}s)")

        if last_projected is not None and projected > last_projected:
            consecutive_increases += 1
            log.write(
                f"    k={k}: projected runtime increased "
                f"({projected:.3f}s > {last_projected:.3f}s); "
                f"consecutive_increases={consecutive_increases}/3\n"
            )
            log.flush()
        else:
            consecutive_increases = 0
        last_projected = projected

        if best_projected is None or projected < best_projected:
            best_projected = projected
            best_k         = k
            best_blocks    = blocks

        if consecutive_increases >= 3:
            msg = (
                f"  [t={t}] Stopping probe sweep after k={k}: "
                "projected runtime increased for 3 consecutive k values."
            )
            print(msg)
            log.write(msg + "\n")
            log.flush()
            break

    # ---- Probe sweep summary ----
    summary_str = {
        k: (f"{v:.3f}s" if v is not None else "pruned/skipped")
        for k, v in probe_summary.items()
    }
    log.write(
        f"\n  t={t}: Probe sweep complete.\n"
        f"    Probe summary (projected totals): {summary_str}\n"
    )

    if best_k is None:
        log.write(f"  t={t}: No valid k found.\n")
        log.flush()
        print(f"  [t={t}] No valid k found.")
        return None, None, None

    log.write(f"    Best k={best_k}, projected={best_projected:.3f}s\n")
    log.flush()
    print(f"\n  [t={t}] Best k={best_k}, projected total={best_projected:.3f}s")

    if probe_only:
        return best_k, best_projected, None

    # ---- Full enumeration on best k ----
    print(f"\n  [t={t}] Running full enumeration on k={best_k} ...")
    full_result, _ = full_run_k(
        best_k, best_blocks, n, all_monomers, mode, all_domains, n_monomers, log
    )

    if full_result is not None and save and save_dir:
        os.makedirs(save_dir, exist_ok=True)
        output_path = os.path.join(
            save_dir, f"hilbert_k{best_k}_t{t}_{mode}.txt"
        )
        save_polymer_vectors(
            full_result["vectors"],
            output_path,
            n_monomers=n_monomers,
            comment=f"covering mode, k={best_k}, t={t}, mode={mode}"
        )

    return best_k, best_projected, full_result


# -------------------------
# Helpers for CLI
# -------------------------

def load_monomers(monomer_file: str) -> list:
    monomers = []
    with open(monomer_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                line = line.split(":", 1)[1].strip()
            if line:
                monomers.append(line)
    return monomers


def validate_args(args, n_monomers, n_domains):
    errors   = []
    warnings = []
    n        = n_monomers if args.mode == "monomer" else n_domains

    if args.t < 1:
        errors.append(f"--t must be >= 1 (got {args.t}).")

    if args.k is None:
        # Probe-sweep mode: need at least one valid k = t+1 .. min(K_MAX, n).
        if args.t >= min(K_MAX, n):
            errors.append(
                f"--t must be < min(K_MAX={K_MAX}, n={n})={min(K_MAX,n)} "
                f"so that at least one valid k = t+1 .. {min(K_MAX,n)} exists."
            )
    else:
        # Fixed-k mode: k must be at least t and at most n.
        if not (args.t <= args.k <= n):
            errors.append(
                f"--k must be in [t={args.t}, n={n}] (got --k {args.k})."
            )

    if not args.fallback_dp and (n >= 100 or args.t > 8):
        issues = []
        if n >= 100:   issues.append(f"n={n} >= 100")
        if args.t > 8: issues.append(f"t={args.t} > 8")
        errors.append(
            f"La Jolla repository limits: n < 100, t <= 8. "
            f"Parameters exceed limits: {', '.join(issues)}. Add --fallback-dp."
        )

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  [!] {w}")
    if errors:
        print("\nErrors:")
        for i, e in enumerate(errors, 1):
            print(f"  [{i}] {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Pareto-Optimal Polymer Enumeration via Hilbert Basis "
                    "(paper Algorithm 1 — covering-design strategy).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Probe-and-prune sweep for t=5, picks the best k automatically\n"
            "  python hilbert_pipeline.py --t 5\n\n"
            "  # Probe-only — estimate per-k runtimes without full enumeration\n"
            "  python hilbert_pipeline.py --t 5 --probe\n\n"
            "  # Fixed k=25, skip the probe sweep\n"
            "  python hilbert_pipeline.py --t 5 --k 25 --save\n\n"
            "  # Domain mode, include full-system base case, save output\n"
            "  python hilbert_pipeline.py --t 4 --mode domain --include-base --save\n\n"
            "Type 's' + Enter during a run to skip the current k."
        )
    )

    parser.add_argument("--t", type=int, default=5, metavar="INT",
                        help="Support bound t (paper Section 5). Default: 5.")
    parser.add_argument("--k", type=int, default=None, metavar="INT",
                        help="If given, skip the probe sweep and run the full "
                             "enumeration at this fixed k. Otherwise the probe "
                             "sweep over k ∈ [t+1, min(25, n)] picks the best k.")
    parser.add_argument("--mode", choices=["monomer", "domain"], default="monomer")
    parser.add_argument("--monomer-file", type=str, default=DEFAULT_MONOMER_FILE,
                        dest="monomer_file")
    parser.add_argument("--include-base", action="store_true", dest="include_base",
                        help="Also probe/run k=n (full system). Ignored when --k is set.")
    parser.add_argument("--fallback-dp", action="store_true", dest="fallback_dp",
                        help="If LJCR has no covering for (n, k, t), build one "
                             "locally via the GPK dynamic-programming "
                             "construction (Gordon-Patashnik-Kuperberg 1995, "
                             "Section 5).")
    parser.add_argument("--probe", action="store_true",
                        help="Probe-only: estimate runtimes without full "
                             "enumeration. Ignored when --k is set.")
    parser.add_argument("--save", action="store_true",
                        help="Save Hilbert basis vectors to disk.")
    parser.add_argument("--save-dir", type=str,
                        default=os.path.join(RESULTS_DIR, "hilbert_output"),
                        dest="save_dir", metavar="PATH")

    args = parser.parse_args()

    if args.save:
        os.makedirs(args.save_dir, exist_ok=True)

    cleanup_normaliz_files()

    all_monomers = load_monomers(args.monomer_file)
    n_monomers   = len(all_monomers)
    all_domains  = get_all_unique_domains(all_monomers)
    n_domains    = len(all_domains)

    print(f"\nLoaded {n_monomers} monomers, {n_domains} unique domain types")
    print(f"Mode: {args.mode}  |  t={args.t}"
          + (f"  |  k={args.k} (fixed)" if args.k is not None else "  |  k=probe-best"))

    validate_args(args, n_monomers, n_domains)

    os.makedirs(LOGS_DIR, exist_ok=True)
    start_input_listener()

    n = n_monomers if args.mode == "monomer" else n_domains

    # ---------- Fixed-k path: skip probe, run full enumeration directly ----
    if args.k is not None:
        log_file = os.path.join(
            LOGS_DIR,
            f"log_fixedk_{args.mode}_n{n_monomers}_t{args.t}_k{args.k}.txt",
        )
        with open(log_file, "a") as log:
            log.write(
                f"Covering (fixed k) — Pareto-Optimal Polymer Enumeration\n"
                f"Started: {datetime.now()}\n"
                f"Mode: {args.mode}  |  t={args.t}  |  k={args.k} (fixed)\n"
                f"n_monomers={n_monomers}  n_domains={n_domains}\n"
                f"fallback_dp={args.fallback_dp}\n"
                + "=" * 70 + "\n"
            )
            try:
                blocks = load_covering_blocks(n, args.k, args.t,
                                              fallback_dp=args.fallback_dp)
            except RuntimeError as e:
                print(f"Could not load covering C({n},{args.k},{args.t}): {e}")
                log.write(f"FAILED to load covering: {e}\n")
                cleanup_normaliz_files()
                return

            print(f"  Loaded {len(blocks)} covering blocks. Running full enumeration ...")
            try:
                full_result, _ = full_run_k(
                    args.k, blocks, n, all_monomers, args.mode,
                    all_domains, n_monomers, log,
                )
            except KeyboardInterrupt:
                print("\nInterrupted."); log.write("\nInterrupted by user.\n")
                log.flush(); cleanup_normaliz_files(); return

            if full_result is not None:
                print(
                    f"\nFixed k={args.k}: wall={full_result['total_wall_time']:.3f}s  "
                    f"normaliz={full_result['total_normaliz_time']:.3f}s  "
                    f"unique_vectors={full_result['unique_vectors']}"
                )
                if args.save:
                    save_polymer_vectors(
                        full_result["vectors"],
                        os.path.join(
                            args.save_dir,
                            f"hilbert_k{args.k}_t{args.t}_{args.mode}.txt",
                        ),
                        n_monomers=n_monomers,
                        comment=f"covering mode (fixed k), k={args.k}, "
                                f"t={args.t}, mode={args.mode}"
                    )
        cleanup_normaliz_files()
        return

    # ---------- Probe-sweep path: run_covering_sweep finds the best k -----
    log_file = os.path.join(
        LOGS_DIR,
        f"log_covering_{args.mode}_n{n_monomers}_t{args.t}"
        f"{'_probe' if args.probe else ''}.txt",
    )
    with open(log_file, "a") as log:
        log.write(
            f"Covering Strategy — Pareto-Optimal Polymer Enumeration\n"
            f"Started: {datetime.now()}\n"
            f"Mode: {args.mode}  |  t={args.t}  |  probe_only={args.probe}\n"
            f"n_monomers={n_monomers}  n_domains={n_domains}\n"
            f"include_base={args.include_base}  fallback_dp={args.fallback_dp}\n"
            + "=" * 70 + "\n"
        )
        try:
            best_k, best_projected, full_result = run_covering_sweep(
                t               = args.t,
                all_monomers    = all_monomers,
                mode            = args.mode,
                log             = log,
                fallback_dp     = args.fallback_dp,
                probe_only      = args.probe,
                include_base    = args.include_base,
                save            = args.save,
                save_dir        = args.save_dir,
            )
        except KeyboardInterrupt:
            print("\nInterrupted."); log.write("\nInterrupted by user.\n")
            log.flush(); cleanup_normaliz_files(); return

    if best_k is not None:
        print(f"\nBest k={best_k}, projected total={best_projected:.3f}s")
        if full_result:
            print(
                f"Full run: wall={full_result['total_wall_time']:.3f}s  "
                f"normaliz={full_result['total_normaliz_time']:.3f}s  "
                f"unique_vectors={full_result['unique_vectors']}"
            )

    cleanup_normaliz_files()

if __name__ == "__main__":
    try:
        main()
    finally:
        cleanup_normaliz_files()
