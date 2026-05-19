from __future__ import annotations

import argparse
import sqlite3
import time
from pathlib import Path
from urllib.request import Request, urlopen

from parse_ljcr_to_sqlite import (
    connect_db as connect_repo_db,
    import_coverdata_csv,
    import_coverdata_json,
    import_covers_csv,
    import_covers_json,
)
from cover_bound_dp_sqlite import (
    connect_db as connect_dp_db,
    derive_single_target,
    derive_upper_bounds_full,
    load_seed_bounds,
    print_target,
    store_run_metadata,
    write_results,
)


def _detect_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".json":
        return "json"
    if suffix == ".csv":
        return "csv"
    raise ValueError(f"Cannot infer file type from {path}. Use .json or .csv files.")


def _count(conn: sqlite3.Connection, table: str) -> int:
    return int(conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])


def _zenodo_url(record: str, filename: str) -> str:
    return f"https://zenodo.org/records/{record}/files/{filename}?download=1"


def _download_file(url: str, dest: Path) -> None:
    """Download `url` to `dest`, resuming a partial .part file when possible."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    existing = part.stat().st_size if part.exists() else 0

    headers = {"User-Agent": "Hilbert-Basis-Algorithm/1.0"}
    if existing:
        headers["Range"] = f"bytes={existing}-"

    req = Request(url, headers=headers)
    mode = "ab" if existing else "wb"
    with urlopen(req, timeout=60) as resp:
        status = getattr(resp, "status", None)
        # Some servers ignore Range and return 200; restart rather than corrupt.
        if existing and status == 200:
            existing = 0
            mode = "wb"
        total_header = resp.headers.get("Content-Length")
        total = int(total_header) + existing if total_header else None
        print(f"      downloading {url}")
        if existing:
            print(f"      resuming at {existing / 1e9:.2f} GB")
        done = existing
        last = time.time()
        with open(part, mode) as f:
            while True:
                chunk = resp.read(8 * 1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                done += len(chunk)
                now = time.time()
                if now - last >= 5:
                    if total:
                        print(f"      {done / 1e9:.2f}/{total / 1e9:.2f} GB ({100 * done / total:.1f}%)")
                    else:
                        print(f"      {done / 1e9:.2f} GB")
                    last = now

    if total is not None and done < total:
        raise RuntimeError(
            f"Download ended early for {dest.name}: got {done} bytes, expected {total}. "
            f"Kept partial file at {part}; rerun to resume."
        )
    part.replace(dest)


def _resolve_ljcr_inputs(args) -> tuple[Path, Path | None]:
    """Return local coverdata and optional covers paths, downloading as needed."""
    if args.coverdata is not None:
        coverdata = args.coverdata
    else:
        download_dir = args.download_dir or (args.out_dir / "downloads")
        coverdata = download_dir / "coverdata.json"
        if not coverdata.exists():
            _download_file(_zenodo_url(args.zenodo_record, "coverdata.json"), coverdata)
        else:
            print(f"      using cached {coverdata}")

    if args.bounds_only:
        return coverdata, None

    if args.covers is not None:
        return coverdata, args.covers

    download_dir = args.download_dir or (args.out_dir / "downloads")
    covers = download_dir / "covers.json"
    if covers.exists():
        # Guard against an interrupted older run that accidentally promoted a
        # partial download to covers.json. A complete top-level JSON object must
        # end with '}' after trailing whitespace.
        try:
            with open(covers, "rb") as f:
                f.seek(max(0, covers.stat().st_size - 4096))
                tail = f.read().rstrip()
            if not tail.endswith(b"}"):
                part = covers.with_suffix(covers.suffix + ".part")
                print(f"      cached {covers} looks incomplete; moving it to {part} for resume")
                if part.exists():
                    part.unlink()
                covers.rename(part)
        except OSError:
            pass

    if not covers.exists():
        _download_file(_zenodo_url(args.zenodo_record, "covers.json"), covers)
    else:
        print(f"      using cached {covers}")

    return coverdata, covers


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build the offline SQLite covering-design databases used by "
            "hilbert_pipeline.py --fallback-dp."
        )
    )
    parser.add_argument("--coverdata", type=Path,
                        help="Local LJCR coverdata file (.json or .csv): sizes/bounds. If omitted, download from Zenodo.")
    parser.add_argument("--covers", type=Path,
                        help="Local LJCR covers file (.json or .csv): actual seed blocks. If omitted, download from Zenodo unless --bounds-only is set.")
    parser.add_argument("--bounds-only", action="store_true",
                        help="Build only the compact DP recipe DB from coverdata sizes. Do not download/import covers.json seed blocks.")
    parser.add_argument("--zenodo-record", default="19735294",
                        help="Zenodo record to download from when --coverdata/--covers are omitted. Default: 19735294.")
    parser.add_argument("--download-dir", type=Path,
                        help="Where downloaded LJCR JSON files are cached. Default: OUT_DIR/downloads.")
    parser.add_argument("--out-dir", type=Path, default=Path("data/covering_design_rebuild"),
                        help="Directory for generated SQLite files.")
    parser.add_argument("--V", type=int, required=True,
                        help="Maximum v to precompute.")
    parser.add_argument("--K", type=int, required=True,
                        help="Maximum k to precompute.")
    parser.add_argument("--T", type=int, default=8,
                        help="Maximum t to precompute. Default: 8.")
    parser.add_argument("--run-name", type=str,
                        help="Name for this DP run. Default: gpk_V{V}_K{K}_T{T}.")
    parser.add_argument("--mode", choices=["full", "target"], default="full",
                        help="Precompute every state up to V/K/T, or only one target state.")
    parser.add_argument("--target", nargs=3, type=int, metavar=("v", "k", "t"),
                        help="Target state for --mode target.")
    parser.add_argument("--replace", action="store_true",
                        help="Overwrite existing imported rows / DP run with the same name.")
    args = parser.parse_args()

    if args.mode == "target" and args.target is None:
        parser.error("--mode target requires --target v k t")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    repo_db = args.out_dir / "ljcr.sqlite"
    dp_db = args.out_dir / "gpk_dp.sqlite"
    run_name = args.run_name or f"gpk_V{args.V}_K{args.K}_T{args.T}"

    coverdata_path, covers_path = _resolve_ljcr_inputs(args)
    coverdata_kind = _detect_kind(coverdata_path)
    covers_kind = _detect_kind(covers_path) if covers_path is not None else None

    print(f"[1/2] Importing LJCR bounds into {repo_db}")
    repo_conn = connect_repo_db(repo_db)
    repo_conn.row_factory = sqlite3.Row
    try:
        if coverdata_kind == "json":
            import_coverdata_json(repo_conn, coverdata_path, replace=args.replace)
        else:
            import_coverdata_csv(repo_conn, coverdata_path, replace=args.replace)

        if covers_path is not None:
            if covers_kind == "json":
                import_covers_json(repo_conn, covers_path, replace=args.replace)
            else:
                import_covers_csv(repo_conn, covers_path, replace=args.replace)

        n_bounds = _count(repo_conn, "seed_bounds")
        n_designs = _count(repo_conn, "seed_design_meta")
        n_blocks = _count(repo_conn, "seed_design_blocks")
        print(f"      seed bounds : {n_bounds}")
        print(f"      seed designs: {n_designs}" + (" (skipped by --bounds-only)" if covers_path is None else ""))
        print(f"      seed blocks : {n_blocks}" + (" (skipped by --bounds-only)" if covers_path is None else ""))

        print(f"[2/2] Computing GPK DP reason table into {dp_db}")
        dp_conn = connect_dp_db(dp_db)
        try:
            target = tuple(args.target) if args.target is not None else None
            store_run_metadata(
                dp_conn,
                run_name,
                args.mode,
                args.V,
                args.K,
                args.T,
                target,
                str(repo_db),
                replace=args.replace,
            )
            seed_bounds = load_seed_bounds(repo_conn)
            if args.mode == "full":
                upper, reasons = derive_upper_bounds_full(args.V, args.K, args.T, seed_bounds)
            else:
                upper, reasons = derive_single_target(target, args.V, args.K, args.T, seed_bounds)  # type: ignore[arg-type]
            write_results(dp_conn, run_name, upper, reasons)
            if target is not None:
                print_target(dp_conn, run_name, target)  # type: ignore[arg-type]
            else:
                print(f"      stored states: {len(reasons)}")
        finally:
            dp_conn.close()
    finally:
        repo_conn.close()

    print("\nDone. To use this database with hilbert_pipeline.py:")
    print(f"  export COVERING_DP_DB={dp_db.resolve()}")
    print(f"  export COVERING_DP_RUN_NAME={run_name}")
    if not args.bounds_only:
        print(f"  export COVERING_REPO_DB={repo_db.resolve()}")
    else:
        print("  # COVERING_REPO_DB is optional in bounds-only mode; seed designs will be fetched online.")
    print("  python src/hilbert_pipeline.py ... --fallback-dp")


if __name__ == "__main__":
    main()
