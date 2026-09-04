#!/usr/bin/env python3
"""
Author: Anna Korda

Cross-conformation compound-ID overlap, over apply_score_cutoff.py's already-built
final_hits_<conf>_cutoff<C>.csv files -- just compound_id set intersections, no
structures/fingerprints needed, so this is fast (seconds) regardless of set size.

Writes a NON-symmetric percentage matrix: matrix[from][to] = what fraction of
`from`'s compounds are also in `to`, i.e. divided by |from|'s own count -- since
c1/ref1/c5 don't have the same number of hits, matrix[from][to] != matrix[to][from]
in general (Anna's own framing: "the overlap we did was a non symmetric matrix...
since they do not have the same number of counts before").

Usage:
    python hit_overlap.py --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/hit_overlap \\
        --conformations c1,ref1,c5 --cutoff -7.0
"""
import argparse
import csv
import time
from pathlib import Path


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def load_id_set(vs_results_dir, conf, cutoff):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found -- run apply_score_cutoff.py for {conf} first")
    with open(csv_path) as f:
        return {row["compound_id"] for row in csv.DictReader(f)}


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

    sets = {}
    for conf in confs:
        sets[conf] = load_id_set(args.vs_results_dir, conf, args.cutoff)
        log(f"{conf}: {len(sets[conf]):,} compound IDs loaded")

    long_rows = []
    for a in confs:
        for b in confs:
            inter = sets[a] if a == b else (sets[a] & sets[b])
            n_a, n_b, n_int = len(sets[a]), len(sets[b]), len(inter)
            pct_of_a_in_b = 100.0 * n_int / n_a if n_a else 0.0
            long_rows.append({
                "from": a, "to": b, "n_from": n_a, "n_to": n_b,
                "n_intersection": n_int, "pct_of_from_in_to": f"{pct_of_a_in_b:.4f}",
            })
            log(f"  {a} -> {b}: {n_int:,} shared ({pct_of_a_in_b:.2f}% of {a}'s {n_a:,})")

    long_path = out_dir / "overlap_counts.csv"
    with open(long_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["from", "to", "n_from", "n_to", "n_intersection", "pct_of_from_in_to"])
        w.writeheader()
        w.writerows(long_rows)
    log(f"Wrote {long_path}")

    matrix_path = out_dir / "overlap_matrix.csv"
    with open(matrix_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + confs)
        for a in confs:
            row = [a]
            for b in confs:
                n_a, n_b = len(sets[a]), len(sets[b])
                n_int = len(sets[a] if a == b else (sets[a] & sets[b]))
                row.append(f"{100.0 * n_int / n_a:.4f}" if n_a else "0.0000")
            w.writerow(row)
    log(f"Wrote {matrix_path} (row-normalized: value = %% of ROW conformation's hits also in COLUMN conformation)")

    summary_path = out_dir / "run_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"Cross-conformation hit overlap -- cutoff {args.cutoff:g}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for conf in confs:
            f.write(f"{conf}: {len(sets[conf]):,} compounds\n")
        f.write("\noverlap_matrix.csv is NON-symmetric: matrix[row][col] = "
                "100 * |row set intersect col set| / |row set|\n")
    log(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
