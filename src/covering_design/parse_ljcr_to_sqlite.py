from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

NAME_RE = re.compile(r"C\((\d+),(\d+),(\d+)\)")


def parse_name(name: str) -> Tuple[int, int, int]:
    m = NAME_RE.fullmatch(name.strip())
    if not m:
        raise ValueError(f"Bad covering name: {name!r}")
    return tuple(map(int, m.groups()))  # type: ignore[return-value]


def connect_db(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA foreign_keys=ON")
    create_schema(conn)
    return conn


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS seed_bounds (
            v INTEGER NOT NULL,
            k INTEGER NOT NULL,
            t INTEGER NOT NULL,
            name TEXT NOT NULL,
            size INTEGER NOT NULL,
            low_bd INTEGER,
            best_method TEXT,
            best_creator TEXT,
            best_timestamp TEXT,
            improvement_count INTEGER,
            PRIMARY KEY (v, k, t)
        );

        CREATE TABLE IF NOT EXISTS seed_design_blocks (
            v INTEGER NOT NULL,
            k INTEGER NOT NULL,
            t INTEGER NOT NULL,
            block_index INTEGER NOT NULL,
            block_size INTEGER NOT NULL,
            block TEXT NOT NULL,
            PRIMARY KEY (v, k, t, block_index)
        );

        CREATE TABLE IF NOT EXISTS seed_design_meta (
            v INTEGER NOT NULL,
            k INTEGER NOT NULL,
            t INTEGER NOT NULL,
            block_count INTEGER NOT NULL,
            PRIMARY KEY (v, k, t)
        );

        CREATE TABLE IF NOT EXISTS import_log (
            import_id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_type TEXT NOT NULL,
            source_path TEXT NOT NULL,
            imported_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
            detail TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_seed_design_blocks_key
        ON seed_design_blocks(v, k, t);
        """
    )
    conn.commit()


def parse_coverdata_record(name: str, rec: dict) -> Tuple[int, int, int, str, int, Optional[int], Optional[str], Optional[str], Optional[str], int]:
    v, k, t = parse_name(name)
    imps = rec.get("imps", []) or []
    best_imp = imps[0] if imps else [None, None, None, None]
    low_bd = rec.get("low_bd")
    return (
        v,
        k,
        t,
        name,
        int(rec["size"]),
        None if low_bd in (None, "") else int(low_bd),
        best_imp[1],
        best_imp[2],
        best_imp[3],
        len(imps),
    )


def import_coverdata_json(conn: sqlite3.Connection, coverdata_json: Path, replace: bool) -> None:
    with open(coverdata_json, "r", encoding="utf-8") as f:
        data: Dict[str, dict] = json.load(f)

    sql = (
        "INSERT OR REPLACE INTO seed_bounds "
        "(v,k,t,name,size,low_bd,best_method,best_creator,best_timestamp,improvement_count) "
        "VALUES (?,?,?,?,?,?,?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_bounds "
        "(v,k,t,name,size,low_bd,best_method,best_creator,best_timestamp,improvement_count) "
        "VALUES (?,?,?,?,?,?,?,?,?,?)"
    )
    rows = [parse_coverdata_record(name, rec) for name, rec in data.items()]
    with conn:
        conn.executemany(sql, rows)
        conn.execute(
            "INSERT INTO import_log(source_type, source_path, detail) VALUES (?,?,?)",
            ("coverdata_json", str(coverdata_json), f"rows={len(rows)}"),
        )


# Streaming parser for huge top-level JSON objects

def iter_top_level_json_object(path: str | Path, chunk_size: int = 1 << 20) -> Iterator[Tuple[str, object]]:
    decoder = json.JSONDecoder()

    with open(path, "r", encoding="utf-8") as f:
        buf = ""
        eof = False

        def fill() -> None:
            nonlocal buf, eof
            if eof:
                return
            chunk = f.read(chunk_size)
            if chunk == "":
                eof = True
            else:
                buf += chunk

        fill()
        if buf.startswith("\ufeff"):
            buf = buf[1:]
        buf = buf.lstrip()
        if not buf:
            raise ValueError("Empty file")
        if not buf.startswith("{"):
            preview = buf[:120].replace("\n", "\\n")
            raise ValueError(f"Expected top-level JSON object. First bytes: {preview!r}")
        buf = buf[1:]

        while True:
            while True:
                buf = buf.lstrip()
                if buf:
                    break
                if eof:
                    raise ValueError("Unexpected EOF while parsing top-level JSON object")
                fill()

            if buf.startswith("}"):
                return

            while True:
                try:
                    key, idx = decoder.raw_decode(buf)
                    break
                except json.JSONDecodeError as e:
                    if eof:
                        preview = buf[:200].replace("\n", "\\n")
                        raise ValueError(
                            "Failed while parsing a top-level key. "
                            "The file is likely truncated or not the expected JSON. "
                            f"Remaining buffer starts with: {preview!r}"
                        ) from e
                    fill()
            if not isinstance(key, str):
                raise ValueError("Expected string key in top-level object")
            buf = buf[idx:].lstrip()

            while not buf.startswith(":"):
                if eof:
                    raise ValueError("Expected ':' after key")
                fill()
                buf = buf.lstrip()
            buf = buf[1:].lstrip()

            while True:
                try:
                    value, idx = decoder.raw_decode(buf)
                    break
                except json.JSONDecodeError as e:
                    if eof:
                        preview = buf[:200].replace("\n", "\\n")
                        raise ValueError(
                            "Failed while parsing a top-level value. "
                            "The file is likely truncated or not the expected JSON. "
                            f"Remaining buffer starts with: {preview!r}"
                        ) from e
                    fill()
            buf = buf[idx:]
            yield key, value

            while True:
                buf = buf.lstrip()
                if buf:
                    break
                if eof:
                    raise ValueError("Unexpected EOF after JSON value")
                fill()

            if buf.startswith(","):
                buf = buf[1:]
                continue
            if buf.startswith("}"):
                return
            if eof:
                preview = buf[:200].replace("\n", "\\n")
                raise ValueError(f"Unexpected trailing content near: {preview!r}")
            fill()


def _extract_blocks(value: object) -> List[List[int]]:
    if not isinstance(value, list):
        raise ValueError(f"Expected list of blocks, got {type(value)!r}")
    blocks: List[List[int]] = []
    for item in value:
        if isinstance(item, list) and all(isinstance(x, int) for x in item):
            blocks.append([int(x) for x in item])
    return blocks


def import_covers_json(conn: sqlite3.Connection, covers_json: Path, replace: bool, batch_size: int = 10000) -> None:
    if replace:
        with conn:
            conn.execute("DELETE FROM seed_design_blocks")
            conn.execute("DELETE FROM seed_design_meta")

    insert_sql = (
        "INSERT OR REPLACE INTO seed_design_blocks (v,k,t,block_index,block_size,block) VALUES (?,?,?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_design_blocks (v,k,t,block_index,block_size,block) VALUES (?,?,?,?,?,?)"
    )
    meta_sql = (
        "INSERT OR REPLACE INTO seed_design_meta (v,k,t,block_count) VALUES (?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_design_meta (v,k,t,block_count) VALUES (?,?,?,?)"
    )

    block_rows: List[Tuple[int, int, int, int, int, str]] = []
    meta_rows: List[Tuple[int, int, int, int]] = []
    count_keys = 0
    count_blocks = 0

    def flush() -> None:
        nonlocal block_rows, meta_rows
        if not block_rows and not meta_rows:
            return
        with conn:
            if block_rows:
                conn.executemany(insert_sql, block_rows)
            if meta_rows:
                conn.executemany(meta_sql, meta_rows)
        block_rows = []
        meta_rows = []

    for name, value in iter_top_level_json_object(covers_json):
        v, k, t = parse_name(name)
        blocks = _extract_blocks(value)
        meta_rows.append((v, k, t, len(blocks)))
        count_keys += 1
        for idx, block in enumerate(blocks):
            block_rows.append((v, k, t, idx, len(block), " ".join(map(str, block))))
            count_blocks += 1
        if len(block_rows) >= batch_size:
            flush()
    flush()
    with conn:
        conn.execute(
            "INSERT INTO import_log(source_type, source_path, detail) VALUES (?,?,?)",
            ("covers_json", str(covers_json), f"keys={count_keys}, blocks={count_blocks}"),
        )


def import_coverdata_csv(conn: sqlite3.Connection, coverdata_csv: Path, replace: bool) -> None:
    sql = (
        "INSERT OR REPLACE INTO seed_bounds (v,k,t,name,size,low_bd,best_method,best_creator,best_timestamp,improvement_count) VALUES (?,?,?,?,?,?,?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_bounds (v,k,t,name,size,low_bd,best_method,best_creator,best_timestamp,improvement_count) VALUES (?,?,?,?,?,?,?,?,?,?)"
    )
    rows = []
    with open(coverdata_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                (
                    int(row["v"]),
                    int(row["k"]),
                    int(row["t"]),
                    row.get("name") or f"C({row['v']},{row['k']},{row['t']})",
                    int(row["size"]),
                    int(row["low_bd"]) if row.get("low_bd") not in (None, "") else None,
                    row.get("best_method"),
                    row.get("best_creator"),
                    row.get("best_timestamp"),
                    int(row.get("improvement_count") or 0),
                )
            )
    with conn:
        conn.executemany(sql, rows)
        conn.execute(
            "INSERT INTO import_log(source_type, source_path, detail) VALUES (?,?,?)",
            ("coverdata_csv", str(coverdata_csv), f"rows={len(rows)}"),
        )


def import_covers_csv(conn: sqlite3.Connection, covers_csv: Path, replace: bool, batch_size: int = 10000) -> None:
    if replace:
        with conn:
            conn.execute("DELETE FROM seed_design_blocks")
            conn.execute("DELETE FROM seed_design_meta")

    insert_sql = (
        "INSERT OR REPLACE INTO seed_design_blocks (v,k,t,block_index,block_size,block) VALUES (?,?,?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_design_blocks (v,k,t,block_index,block_size,block) VALUES (?,?,?,?,?,?)"
    )

    block_rows: List[Tuple[int, int, int, int, int, str]] = []
    meta_count: Dict[Tuple[int, int, int], int] = {}

    def flush() -> None:
        nonlocal block_rows
        if not block_rows:
            return
        with conn:
            conn.executemany(insert_sql, block_rows)
        block_rows = []

    with open(covers_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (int(row["v"]), int(row["k"]), int(row["t"]))
            idx = int(row.get("block_index") or meta_count.get(key, 0))
            block = row["block"].strip()
            block_size = int(row.get("block_size") or len(block.split()))
            block_rows.append((key[0], key[1], key[2], idx, block_size, block))
            meta_count[key] = meta_count.get(key, 0) + 1
            if len(block_rows) >= batch_size:
                flush()
    flush()

    meta_sql = (
        "INSERT OR REPLACE INTO seed_design_meta (v,k,t,block_count) VALUES (?,?,?,?)"
        if replace
        else
        "INSERT OR IGNORE INTO seed_design_meta (v,k,t,block_count) VALUES (?,?,?,?)"
    )
    with conn:
        conn.executemany(meta_sql, [(v, k, t, c) for (v, k, t), c in meta_count.items()])
        conn.execute(
            "INSERT INTO import_log(source_type, source_path, detail) VALUES (?,?,?)",
            ("covers_csv", str(covers_csv), f"keys={len(meta_count)}, blocks={sum(meta_count.values())}"),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Import La Jolla Coverings Repository data into an indexed SQLite database.")
    parser.add_argument("--db", type=Path, required=True, help="Output SQLite database path")
    parser.add_argument("--coverdata-json", type=Path)
    parser.add_argument("--coverdata-csv", type=Path)
    parser.add_argument("--covers-json", type=Path)
    parser.add_argument("--covers-csv", type=Path)
    parser.add_argument("--replace", action="store_true", help="Overwrite existing rows for imported tables")
    args = parser.parse_args()

    if not any([args.coverdata_json, args.coverdata_csv, args.covers_json, args.covers_csv]):
        parser.error("Provide at least one of --coverdata-json, --coverdata-csv, --covers-json, --covers-csv")

    conn = connect_db(args.db)
    try:
        if args.coverdata_json:
            import_coverdata_json(conn, args.coverdata_json, replace=args.replace)
            print(f"Imported numeric bounds from {args.coverdata_json} into {args.db}")
        if args.coverdata_csv:
            import_coverdata_csv(conn, args.coverdata_csv, replace=args.replace)
            print(f"Imported numeric bounds from {args.coverdata_csv} into {args.db}")
        if args.covers_json:
            import_covers_json(conn, args.covers_json, replace=args.replace)
            print(f"Imported covering designs from {args.covers_json} into {args.db}")
        if args.covers_csv:
            import_covers_csv(conn, args.covers_csv, replace=args.replace)
            print(f"Imported covering designs from {args.covers_csv} into {args.db}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
