#!/usr/bin/env python3
"""
Author: Anna Korda

Takes collect_final_hits.py's already-built, already-rank-ordered output
and keeps only compounds with ad4_score <= --cutoff (default -6.0 kcal/mol
-- some compounds pass the salt-bridge + strain filters with a positive AD4
score, i.e. the scoring function itself predicts unfavorable binding; this
cuts those out).

Does NOT re-fetch anything from the original per-compound files. Since
final_hits_<conf>.csv is already sorted best-score-first and the
final_hits_<conf>_partNN.sdf files are in that same order, the cutoff
selects a contiguous prefix -- this just reads the already-built part
files (fast, local, no GPFS scatter) and re-splits that prefix into new
output files. Molecules are never touched, only reorganized -- same
guarantee as every other script in this pipeline.

Usage:
    python apply_score_cutoff.py --conformation c1 \\
        --vs_results_dir vs_results \\
        --cutoff -6.0 --chunk_size 50000
"""
import argparse
import csv
from pathlib import Path


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-6.0,
                     help="keep compounds with ad4_score <= this value")
    ap.add_argument("--chunk_size", type=int, default=50000)
    args = ap.parse_args()

    conf = args.conformation
    src_dir = Path(args.vs_results_dir) / f"final_hits_{conf}"
    csv_path = src_dir / f"final_hits_{conf}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found -- run collect_final_hits.py for this conformation first")

    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    n_total = len(rows)

    # rows are already sorted best-score-first (collect_final_hits.py's own
    # guarantee) -- verify that assumption rather than trust it blindly,
    # since everything below depends on it.
    scores = [float(r["ad4_score"]) for r in rows]
    if scores != sorted(scores):
        raise SystemExit(f"{csv_path} is not sorted by ad4_score ascending as expected -- "
                         f"refusing to guess a prefix boundary. Check the file.")

    n_keep = sum(1 for s in scores if s <= args.cutoff)
    print(f"{n_total:,} total compounds in {csv_path.name}; "
          f"{n_keep:,} have ad4_score <= {args.cutoff} (kept, a contiguous best-scoring prefix)")
    if n_keep == 0:
        raise SystemExit(f"No compounds meet the cutoff of {args.cutoff} -- nothing to write.")

    out_dir = Path(args.vs_results_dir) / f"final_hits_{conf}_cutoff{args.cutoff:g}"
    out_dir.mkdir(parents=True, exist_ok=True)

    kept_rows = rows[:n_keep]
    with open(out_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "compound_id", "ad4_score", "total_TEU", "single_TEU"])
        w.writeheader()
        w.writerows(kept_rows)

    part_files = sorted(src_dir.glob(f"final_hits_{conf}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No final_hits_{conf}_part*.sdf files found in {src_dir}")

    part_num = 1
    count_in_part = 0
    written = 0
    out_f = open(out_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}_part{part_num:02d}.sdf", "w")
    for pf in part_files:
        if written >= n_keep:
            break
        for block in _split_sdf_blocks(pf.read_text()):
            if written >= n_keep:
                break
            if count_in_part >= args.chunk_size:
                out_f.close()
                part_num += 1
                count_in_part = 0
                out_f = open(out_dir / f"final_hits_{conf}_cutoff{args.cutoff:g}_part{part_num:02d}.sdf", "w")
            out_f.write(block)
            count_in_part += 1
            written += 1
    out_f.close()

    if written != n_keep:
        raise SystemExit(f"Wrote {written:,} but expected {n_keep:,} -- the part files may not actually "
                         f"be in the same order as the CSV, don't trust this output, investigate.")

    print(f"Wrote {written:,} compounds across {part_num} file(s) to {out_dir}/")


if __name__ == "__main__":
    main()
