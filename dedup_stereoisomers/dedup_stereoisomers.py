#!/usr/bin/env python3
"""
Author: Anna Korda

Deduplicates apply_score_cutoff.py's final_hits_<conf>_cutoff<C>.csv by base
compound: Enamine REAL Space compound_ids that share everything before a
trailing `_S<n>_T<m>` suffix are the same 2D compound -- different
stereoisomer/tautomer/protonation-state enumerations (confirmed real, from
the `cxcalc`-generated REMARK block in these SDFs), not different compounds.
Real, verified consequence: with this project's default ECFP4 fingerprint
settings (no includeChirality), a real R/S enantiomer pair comes back at
Tanimoto = 1.0 -- completely indistinguishable -- so keeping every
stereoisomer inflates both the raw hit count and Butina cluster membership
with entries that add zero diversity signal.

Keeps exactly one entry per base compound: whichever has the best (lowest,
most negative) AD4 score -- same "best AD4 score wins" convention already
used everywhere else in this pipeline (salt-bridge pose selection,
final_hits ranking). Writes a new CSV (re-ranked by ad4_score ascending) and
the corresponding SDF poses, pure text-level splice -- same bytes as the
original blocks, never touched, never opened with a chemistry library.

Usage:
    python dedup_stereoisomers.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --out_dir /users/gpcr/annak/ultra-large/vs_results/dedup \\
        --cutoff -7.0
"""
import argparse
import csv
import re
import time
from pathlib import Path

BASE_ID_RE = re.compile(r"^(.*)_S\d+_T\d+$")


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def _split_sdf_blocks(text):
    blocks = []
    start = 0
    idx = text.find("$$$$")
    while idx != -1:
        end = idx + 4
        if text[end:end + 2] == "\r\n":
            end += 2
        elif text[end:end + 1] == "\n":
            end += 1
        blocks.append(text[start:end])
        start = end
        idx = text.find("$$$$", start)
    if text[start:].strip():
        blocks.append(text[start:])
    return blocks


def base_id(compound_id):
    m = BASE_ID_RE.match(compound_id)
    return m.group(1) if m else compound_id


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    args = ap.parse_args()

    conf = args.conformation
    hit_dir = Path(args.vs_results_dir) / f"final_hits_{conf}_cutoff{args.cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    log(f"{conf}: {len(rows):,} total compound_ids loaded")

    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{args.cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(rows):
        raise SystemExit(f"{conf}: {len(rows)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order, refusing to guess a correspondence")

    # group by base compound, keep the row with the best (lowest) ad4_score
    best_by_base = {}
    block_by_id = {}
    for row, block in zip(rows, blocks):
        cid = row["compound_id"]
        block_by_id[cid] = block
        base = base_id(cid)
        score = float(row["ad4_score"])
        if base not in best_by_base or score < float(best_by_base[base]["ad4_score"]):
            best_by_base[base] = row

    kept_rows = sorted(best_by_base.values(), key=lambda r: float(r["ad4_score"]))
    n_total, n_kept = len(rows), len(kept_rows)
    n_redundant = n_total - n_kept
    log(f"{conf}: {n_total:,} total -> {n_kept:,} unique base compounds "
        f"({n_redundant:,} redundant stereo/tautomer entries removed, "
        f"{100 * n_redundant / n_total:.1f}%)")

    out_dir = Path(args.out_dir) / f"final_hits_{conf}_cutoff{args.cutoff:g}_dedup"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path_out = out_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}_dedup.csv"
    with open(csv_path_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "compound_id", "ad4_score", "total_TEU", "single_TEU"])
        w.writeheader()
        for rank, row in enumerate(kept_rows, start=1):
            w.writerow({"rank": rank, "compound_id": row["compound_id"], "ad4_score": row["ad4_score"],
                       "total_TEU": row["total_TEU"], "single_TEU": row["single_TEU"]})
    log(f"Wrote {csv_path_out}")

    sdf_path_out = out_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}_dedup.sdf"
    with open(sdf_path_out, "w") as f:
        for row in kept_rows:
            f.write(block_by_id[row["compound_id"]])
    log(f"Wrote {sdf_path_out}")

    return n_total, n_kept, n_redundant


if __name__ == "__main__":
    main()
