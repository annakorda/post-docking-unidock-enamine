#!/usr/bin/env python3
"""
Author: Anna Korda

Compounds present in ALL requested conformations' final hit sets (default
c1+ref1+c5) -- the actual intersection this time, not just the pairwise
percentages hit_overlap.py reports. Extracts each overlapping compound's own
winning-pose block from every conformation's already-built
final_hits_<conf>_cutoff<C>_partNN.sdf (byte-identical, same "never touch a
docked molecule" guarantee as the rest of this repo -- pure text-level lookup
by compound_id, blocks are never opened with a chemistry library here), and
reports the real count.

Writes one SDF per conformation (triple_overlap_<conf>.sdf) so the same
compound's pose in each conformation can be loaded side by side in Maestro --
row N in one file is the same compound as row N in the others, ordered by
combined AD4 score (sum across conformations, best first).

Usage:
    python triple_overlap.py --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/triple_overlap \\
        --conformations c1,ref1,c5 --cutoff -7.0
"""
import argparse
import csv
import time
from functools import reduce
from pathlib import Path


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


def load_hit_data(vs_results_dir, conf, cutoff):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found -- run apply_score_cutoff.py for {conf} first")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    ids = [row["compound_id"] for row in rows]
    scores = {row["compound_id"]: float(row["ad4_score"]) for row in rows}

    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order (apply_score_cutoff.py's own guarantee), "
                          f"refusing to guess a correspondence")

    id_to_block = dict(zip(ids, blocks))
    return set(ids), id_to_block, scores


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--conformations", default="c1,ref1,c5")
    ap.add_argument("--cutoff", type=float, default=-7.0)
    args = ap.parse_args()

    confs = args.conformations.split(",")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    id_sets, id_to_block, scores = {}, {}, {}
    for conf in confs:
        log(f"Loading {conf} hits (cutoff {args.cutoff:g}) ...")
        id_sets[conf], id_to_block[conf], scores[conf] = load_hit_data(args.vs_results_dir, conf, args.cutoff)
        log(f"  {len(id_sets[conf]):,} compounds")

    overlap_ids = reduce(lambda a, b: a & b, id_sets.values())
    log(f"{len(overlap_ids):,} compound(s) present in ALL {len(confs)} conformations "
        f"({'+'.join(confs)})")

    ranked = sorted(overlap_ids, key=lambda cid: sum(scores[conf][cid] for conf in confs))

    for conf in confs:
        out_path = out_dir / f"triple_overlap_{conf}.sdf"
        with open(out_path, "w") as f:
            for cid in ranked:
                f.write(id_to_block[conf][cid])
        log(f"Wrote {out_path} ({len(ranked):,} compounds, same order as the other conformations' files)")

    report_path = out_dir / "triple_overlap_report.csv"
    with open(report_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound_id"] + [f"ad4_score_{conf}" for conf in confs] + ["ad4_score_sum"])
        for rank, cid in enumerate(ranked, start=1):
            row_scores = [scores[conf][cid] for conf in confs]
            w.writerow([rank, cid] + [f"{s:.4f}" for s in row_scores] + [f"{sum(row_scores):.4f}"])
    log(f"Wrote {report_path}")

    summary_path = out_dir / "run_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"Triple overlap -- conformations {'+'.join(confs)}, cutoff {args.cutoff:g}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for conf in confs:
            f.write(f"{conf}: {len(id_sets[conf]):,} compounds\n")
        f.write(f"\nIn ALL {len(confs)} conformations: {len(overlap_ids):,} compounds\n")
        f.write(f"\ntriple_overlap_<conf>.sdf files are row-aligned (same compound, same row, "
                f"across every conformation's file) and ranked by summed AD4 score (best first) -- "
                f"see triple_overlap_report.csv for the per-conformation scores.\n")
    log(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
