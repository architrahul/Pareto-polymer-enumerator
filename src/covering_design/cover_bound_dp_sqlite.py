from __future__ import annotations

import argparse
import json
import math
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

Key = Tuple[int, int, int]
INF = 10**30


@dataclass(frozen=True)
class SeedReason:
    source: str
    detail: str


@dataclass(frozen=True)
class LiftReason:
    src: Key


@dataclass(frozen=True)
class WeakenReason:
    src: Key


@dataclass(frozen=True)
class IntervalDirectNode:
    i: int
    j: int
    ell: int
    left: Key
    right: Key
    cost: int


@dataclass(frozen=True)
class IntervalSplitNode:
    i: int
    j: int
    r: int
    left: "IntervalNode"
    right: "IntervalNode"
    cost: int


IntervalNode = Union[IntervalDirectNode, IntervalSplitNode]


@dataclass(frozen=True)
class Section5Reason:
    v1: int
    v2: int
    tree: IntervalNode


Reason = Union[SeedReason, LiftReason, WeakenReason, Section5Reason]


def connect_db(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.row_factory = sqlite3.Row
    create_schema(conn)
    return conn


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS dp_runs (
            run_name TEXT PRIMARY KEY,
            mode TEXT NOT NULL,
            V INTEGER NOT NULL,
            K INTEGER NOT NULL,
            T INTEGER NOT NULL,
            target_v INTEGER,
            target_k INTEGER,
            target_t INTEGER,
            repo_db_path TEXT,
            created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS dp_bounds (
            run_name TEXT NOT NULL,
            v INTEGER NOT NULL,
            k INTEGER NOT NULL,
            t INTEGER NOT NULL,
            bound INTEGER NOT NULL,
            reason_kind TEXT NOT NULL,
            source TEXT,
            detail TEXT,
            src_v INTEGER,
            src_k INTEGER,
            src_t INTEGER,
            split_v1 INTEGER,
            split_v2 INTEGER,
            tree_json TEXT,
            PRIMARY KEY (run_name, v, k, t),
            FOREIGN KEY (run_name) REFERENCES dp_runs(run_name) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_dp_bounds_lookup
        ON dp_bounds(run_name, v, k, t);
        """
    )
    conn.commit()


def load_seed_bounds(repo_conn: sqlite3.Connection) -> Dict[Key, int]:
    out: Dict[Key, int] = {}
    for row in repo_conn.execute("SELECT v,k,t,size FROM seed_bounds"):
        out[(row["v"], row["k"], row["t"])] = row["size"]
    return out


def trivial_upper_bound(v: int, k: int, t: int) -> Tuple[int, Optional[str]]:
    if not (0 <= t <= k <= v):
        return INF, None
    if t == 0:
        return 1, "t=0"
    if k == v:
        return 1, "k=v"
    if t == 1:
        return (v + k - 1) // k, "t=1"
    if t == k:
        return math.comb(v, k), "t=k"
    return INF, None


def _build_interval_tree(choice: List[List[object]], i: int, j: int) -> IntervalNode:
    tag, *payload = choice[i][j]
    if tag == "direct":
        ell, left_key, right_key, cost = payload
        return IntervalDirectNode(i=i, j=j, ell=ell, left=left_key, right=right_key, cost=cost)
    r, cost = payload
    return IntervalSplitNode(
        i=i,
        j=j,
        r=r,
        left=_build_interval_tree(choice, i, r),
        right=_build_interval_tree(choice, r + 1, j),
        cost=cost,
    )


def best_section5_for_state(v: int, k: int, t: int, upper: Dict[Key, int]) -> Tuple[int, Optional[Section5Reason]]:
    best = INF
    best_reason: Optional[Section5Reason] = None

    for v1 in range(1, v // 2 + 1):
        v2 = v - v1
        c = [[INF] * (t + 1) for _ in range(t + 1)]
        choice: List[List[object]] = [[None] * (t + 1) for _ in range(t + 1)]

        for length in range(1, t + 2):
            for i in range(0, t - length + 2):
                j = i + length - 1
                state_best = INF
                state_choice = None

                ell_lo = max(j, k - v2, 0)
                ell_hi = min(v1, k - (t - i), k)
                for ell in range(ell_lo, ell_hi + 1):
                    left_key = (v1, ell, j)
                    right_key = (v2, k - ell, t - i)
                    a = upper.get(left_key, INF)
                    b = upper.get(right_key, INF)
                    if a >= INF or b >= INF:
                        continue
                    cand = a * b
                    if cand < state_best:
                        state_best = cand
                        state_choice = ("direct", ell, left_key, right_key, cand)

                for r in range(i, j):
                    a = c[i][r]
                    b = c[r + 1][j]
                    if a >= INF or b >= INF:
                        continue
                    cand = a + b
                    if cand < state_best:
                        state_best = cand
                        state_choice = ("split", r, cand)

                c[i][j] = state_best
                choice[i][j] = state_choice

        if c[0][t] < best:
            best = c[0][t]
            best_reason = Section5Reason(v1=v1, v2=v2, tree=_build_interval_tree(choice, 0, t))

    return best, best_reason


def _node_to_dict(node: IntervalNode) -> dict:
    if isinstance(node, IntervalDirectNode):
        return {
            "type": "direct",
            "i": node.i,
            "j": node.j,
            "ell": node.ell,
            "left": list(node.left),
            "right": list(node.right),
            "cost": node.cost,
        }
    return {
        "type": "split",
        "i": node.i,
        "j": node.j,
        "r": node.r,
        "cost": node.cost,
        "left": _node_to_dict(node.left),
        "right": _node_to_dict(node.right),
    }


def reason_to_row(run_name: str, key: Key, bound: int, reason: Reason) -> Tuple:
    v, k, t = key
    if isinstance(reason, SeedReason):
        return (run_name, v, k, t, bound, "seed", reason.source, reason.detail, None, None, None, None, None, None)
    if isinstance(reason, LiftReason):
        return (run_name, v, k, t, bound, "lift", None, None, reason.src[0], reason.src[1], reason.src[2], None, None, None)
    if isinstance(reason, WeakenReason):
        return (run_name, v, k, t, bound, "weaken", None, None, reason.src[0], reason.src[1], reason.src[2], None, None, None)
    if isinstance(reason, Section5Reason):
        return (
            run_name,
            v,
            k,
            t,
            bound,
            "section5",
            None,
            None,
            None,
            None,
            None,
            reason.v1,
            reason.v2,
            json.dumps(_node_to_dict(reason.tree), separators=(",", ":")),
        )
    raise TypeError(reason)


def store_run_metadata(conn: sqlite3.Connection, run_name: str, mode: str, V: int, K: int, T: int, target: Optional[Key], repo_db_path: str, replace: bool) -> None:
    if replace:
        with conn:
            conn.execute("DELETE FROM dp_bounds WHERE run_name = ?", (run_name,))
            conn.execute("DELETE FROM dp_runs WHERE run_name = ?", (run_name,))
    with conn:
        conn.execute(
            "INSERT INTO dp_runs(run_name, mode, V, K, T, target_v, target_k, target_t, repo_db_path) VALUES (?,?,?,?,?,?,?,?,?)",
            (run_name, mode, V, K, T, None if target is None else target[0], None if target is None else target[1], None if target is None else target[2], repo_db_path),
        )


def derive_upper_bounds_full(V: int, K: int, T: int, seed_bounds: Dict[Key, int]) -> Tuple[Dict[Key, int], Dict[Key, Reason]]:
    upper: Dict[Key, int] = {}
    reason: Dict[Key, Reason] = {}

    for v in range(V + 1):
        for k in range(min(K, v) + 1):
            for t in range(min(T, k), -1, -1):
                key = (v, k, t)
                best = INF
                why: Optional[Reason] = None

                seed_val = seed_bounds.get(key, INF)
                if seed_val < best:
                    best = seed_val
                    why = SeedReason(source="repository", detail="seed size from repository db")

                triv, detail = trivial_upper_bound(v, k, t)
                if triv < best:
                    best = triv
                    why = SeedReason(source="trivial", detail=detail or "")

                if t + 1 <= min(T, k):
                    src = (v, k, t + 1)
                    val = upper[src]
                    if val < best:
                        best = val
                        why = WeakenReason(src=src)

                if v >= 1 and k >= 1:
                    src = (v - 1, k - 1, t)
                    val = upper.get(src, INF)
                    if val < best:
                        best = val
                        why = LiftReason(src=src)

                if t >= 1 and v >= 2 and k >= 1:
                    cand, sec_reason = best_section5_for_state(v, k, t, upper)
                    if cand < best and sec_reason is not None:
                        best = cand
                        why = sec_reason

                upper[key] = best
                if why is None:
                    why = SeedReason(source="unknown", detail="no rule found")
                reason[key] = why

    return upper, reason


def derive_single_target(target: Key, V: int, K: int, T: int, seed_bounds: Dict[Key, int]) -> Tuple[Dict[Key, int], Dict[Key, Reason]]:
    memo_val: Dict[Key, int] = {}
    memo_reason: Dict[Key, Reason] = {}

    def solve(key: Key) -> int:
        if key in memo_val:
            return memo_val[key]
        v, k, t = key
        if not (0 <= v <= V and 0 <= k <= K and 0 <= t <= T and t <= k <= v):
            return INF

        best = INF
        why: Optional[Reason] = None

        seed_val = seed_bounds.get(key, INF)
        if seed_val < best:
            best = seed_val
            why = SeedReason(source="repository", detail="seed size from repository db")

        triv, detail = trivial_upper_bound(v, k, t)
        if triv < best:
            best = triv
            why = SeedReason(source="trivial", detail=detail or "")

        if t + 1 <= min(T, k):
            src = (v, k, t + 1)
            val = solve(src)
            if val < best:
                best = val
                why = WeakenReason(src=src)

        if v >= 1 and k >= 1:
            src = (v - 1, k - 1, t)
            val = solve(src)
            if val < best:
                best = val
                why = LiftReason(src=src)

        if t >= 1 and v >= 2 and k >= 1:
            upper_view = MemoUpperView(solve)
            cand, sec_reason = best_section5_for_state(v, k, t, upper_view)  # type: ignore[arg-type]
            if cand < best and sec_reason is not None:
                best = cand
                why = sec_reason

        memo_val[key] = best
        if why is None:
            why = SeedReason(source="unknown", detail="no rule found")
        memo_reason[key] = why
        return best

    solve(target)
    return memo_val, memo_reason


class MemoUpperView(dict):
    def __init__(self, solve_fn):
        self.solve_fn = solve_fn

    def get(self, key, default=None):
        val = self.solve_fn(key)
        if val >= INF:
            return default if default is not None else INF
        return val

    def __getitem__(self, key):
        val = self.solve_fn(key)
        if val >= INF:
            raise KeyError(key)
        return val


def write_results(conn: sqlite3.Connection, run_name: str, upper: Dict[Key, int], reasons: Dict[Key, Reason]) -> None:
    rows = [reason_to_row(run_name, key, upper[key], reasons[key]) for key in reasons]
    with conn:
        conn.executemany(
            "INSERT OR REPLACE INTO dp_bounds (run_name,v,k,t,bound,reason_kind,source,detail,src_v,src_k,src_t,split_v1,split_v2,tree_json) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            rows,
        )


def print_target(conn: sqlite3.Connection, run_name: str, target: Key) -> None:
    row = conn.execute(
        "SELECT bound, reason_kind, source, detail, src_v, src_k, src_t, split_v1, split_v2 FROM dp_bounds WHERE run_name=? AND v=? AND k=? AND t=?",
        (run_name, target[0], target[1], target[2]),
    ).fetchone()
    if row is None:
        print(f"Target C{target} was not computed in run {run_name!r}")
        return
    print(f"C{target} <= {row['bound']}")
    if row["reason_kind"] == "seed":
        print(f"Reason: {row['source']} ({row['detail']})")
    elif row["reason_kind"] in {"lift", "weaken"}:
        print(f"Reason: {row['reason_kind']} from ({row['src_v']},{row['src_k']},{row['src_t']})")
    elif row["reason_kind"] == "section5":
        print(f"Reason: Section 5 split v={row['split_v1']}+{row['split_v2']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute covering upper bounds and store compact reasons in SQLite.")
    parser.add_argument("--repo-db", type=Path, required=True, help="SQLite DB containing seed_bounds (and optionally seed designs)")
    parser.add_argument("--out-db", type=Path, required=True, help="SQLite DB to store DP runs and bounds")
    parser.add_argument("--run-name", type=str, required=True, help="Logical name for this preprocessing run")
    parser.add_argument("--V", type=int, required=True)
    parser.add_argument("--K", type=int, required=True)
    parser.add_argument("--T", type=int, required=True)
    parser.add_argument("--mode", choices=["full", "target"], default="full")
    parser.add_argument("--target", nargs=3, type=int, metavar=("v", "k", "t"))
    parser.add_argument("--replace", action="store_true", help="Replace an existing run with the same name")
    args = parser.parse_args()

    target = tuple(args.target) if args.target is not None else None
    if args.mode == "target" and target is None:
        parser.error("--mode target requires --target v k t")
    if target is not None:
        tv, tk, tt = target
        if not (0 <= tv <= args.V and 0 <= tk <= args.K and 0 <= tt <= args.T and tt <= tk <= tv):
            parser.error("--target must satisfy 0 <= t <= k <= v and lie within the preprocessing caps V,K,T")

    repo_conn = connect_db(args.repo_db)
    out_conn = connect_db(args.out_db)
    try:
        seed_bounds = load_seed_bounds(repo_conn)
        store_run_metadata(out_conn, args.run_name, args.mode, args.V, args.K, args.T, target, str(args.repo_db), replace=args.replace)
        if args.mode == "full":
            upper, reasons = derive_upper_bounds_full(args.V, args.K, args.T, seed_bounds)
        else:
            upper, reasons = derive_single_target(target, args.V, args.K, args.T, seed_bounds)  # type: ignore[arg-type]
        write_results(out_conn, args.run_name, upper, reasons)
        if target is not None:
            print_target(out_conn, args.run_name, target)
        else:
            print(f"Stored {len(reasons)} states in run {args.run_name!r}")
    finally:
        repo_conn.close()
        out_conn.close()


if __name__ == "__main__":
    main()
