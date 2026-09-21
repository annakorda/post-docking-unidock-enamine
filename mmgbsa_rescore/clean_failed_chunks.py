#!/usr/bin/env python3
"""
Author: Anna Korda

Deletes any mmgbsa_part_*.csv chunk file that contains a real failure (a row
with a non-empty "error" column) -- run before resubmitting mmgbsa_rescore.py
so those compounds get genuinely recomputed instead of resumability blindly
trusting failures/timeouts left over from an earlier broken run (e.g. the
64-worker stall that timed out many compounds under real contention -- those
timeouts are stale, not permanent, once run under working settings).

Chunk files with NO failures are left alone (real, valid completed work,
no need to redo it).

Usage:
    python clean_failed_chunks.py --out_dir results --conformation c1
"""
import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", required=True, help="mmgbsa_rescore.py's --out_dir")
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--dry_run", action="store_true", help="report what would be deleted, don't delete")
    args = ap.parse_args()

    conf_dir = Path(args.out_dir) / f"mmgbsa_{args.conformation}"
    part_files = sorted(conf_dir.glob("mmgbsa_part_*.csv"))
    if not part_files:
        raise SystemExit(f"No mmgbsa_part_*.csv files found in {conf_dir}")

    n_deleted, n_kept, n_failed_rows, n_ok_rows = 0, 0, 0, 0
    for pf in part_files:
        with open(pf) as f:
            rows = list(csv.DictReader(f))
        has_failure = any(r.get("error") for r in rows)
        n_failed_rows += sum(1 for r in rows if r.get("error"))
        n_ok_rows += sum(1 for r in rows if not r.get("error"))
        if has_failure:
            n_deleted += 1
            if not args.dry_run:
                pf.unlink()
        else:
            n_kept += 1

    verb = "Would delete" if args.dry_run else "Deleted"
    print(f"{verb} {n_deleted} chunk file(s) with a real failure inside "
          f"({n_failed_rows:,} failed compounds across them), kept {n_kept} clean chunk file(s) "
          f"({n_ok_rows:,} real successes preserved)")
    if n_deleted and not args.dry_run:
        print("Re-run mmgbsa_rescore.py for this conformation -- the deleted chunks' "
              "compounds will be recomputed, everything else stays resumed.")


if __name__ == "__main__":
    main()
